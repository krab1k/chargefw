#include "cli_support.h"

#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <variant>

namespace chargefw::cli {
namespace {

template <typename Reader>
auto read_collection(Reader& reader, const std::string& input_path,
                     const ImportedExportContext::Format format) -> ImportedCollection {
    std::vector<core::Molecule> molecules;
    std::vector<adapters::ImportedMoleculeRecord> records;

    while (const auto record = reader.next()) {
        auto imported = std::move(*record);
        molecules.push_back(imported.molecule);
        records.push_back(std::move(imported));
    }

    if (molecules.empty()) {
        throw std::runtime_error{"No molecules found in input file: " + input_path};
    }

    return ImportedCollection{.molecules =
                                  core::MoleculeCollection{std::move(molecules), input_path},
                              .export_context = {.records = std::move(records), .format = format}};
}

[[nodiscard]] auto parse_record_selection(const std::string& value)
    -> adapters::gemmi::RecordSelection {
    if (value == "all") {
        return adapters::gemmi::RecordSelection::all;
    }
    if (value == "polymers-and-ligands") {
        return adapters::gemmi::RecordSelection::polymers_and_ligands;
    }
    if (value == "polymers") {
        return adapters::gemmi::RecordSelection::polymers;
    }
    throw std::runtime_error{"Unsupported structural record selection: " + value};
}

[[nodiscard]] auto parse_conformer_selection(const std::string& value)
    -> adapters::ConformerSelection {
    if (value == "first") {
        return adapters::ConformerSelection::first;
    }
    if (value == "all") {
        return adapters::ConformerSelection::all;
    }
    throw std::runtime_error{"Unsupported conformer selection: " + value};
}

[[nodiscard]] auto parse_bond_strategy(const std::string& value) -> adapters::gemmi::BondStrategy {
    if (value == "none") {
        return adapters::gemmi::BondStrategy::none;
    }
    if (value == "explicit") {
        return adapters::gemmi::BondStrategy::explicit_bonds;
    }
    if (value == "templates") {
        return adapters::gemmi::BondStrategy::templates;
    }
    if (value == "hybrid") {
        return adapters::gemmi::BondStrategy::hybrid;
    }
    throw std::runtime_error{"Unsupported structural bond strategy: " + value};
}

[[nodiscard]] auto parse_atom_threshold(const std::string& value, const std::string_view name)
    -> std::optional<std::size_t> {
    if (value == "unlimited") {
        return std::nullopt;
    }

    std::size_t threshold = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), threshold);
    if (error != std::errc{} || end != value.data() + value.size()) {
        throw std::runtime_error{std::string{name} +
                                 " must be a non-negative integer or 'unlimited'"};
    }
    return threshold;
}

[[nodiscard]] auto parse_method_option(const std::string& text)
    -> std::pair<std::string, std::pair<std::string, methods::MethodOptionValue>> {
    const auto equals = text.find('=');
    const auto dot = text.find('.');
    if (dot == std::string::npos || equals == std::string::npos || dot == 0 || dot + 1 >= equals ||
        equals + 1 >= text.size()) {
        throw std::runtime_error{"Method option must have the form METHOD.OPTION=VALUE: " + text};
    }

    const auto method_id = text.substr(0, dot);
    const auto option_id = text.substr(dot + 1, equals - dot - 1);
    const auto value = text.substr(equals + 1);
    const auto* method = methods::method_registry().find(method_id);
    if (method == nullptr) {
        throw std::runtime_error{"Unknown method in method option: " + method_id};
    }
    const auto schema = method->option_schema();
    const auto spec =
        std::ranges::find_if(schema, [&option_id](const methods::MethodOptionSpec& item) {
            return item.id == option_id;
        });
    if (spec == schema.end()) {
        throw std::runtime_error{"Unknown option '" + option_id + "' for method '" + method_id +
                                 "'"};
    }

    switch (spec->type) {
    case methods::MethodOptionType::boolean:
        if (value == "true" || value == "1") {
            return {method_id, {option_id, true}};
        }
        if (value == "false" || value == "0") {
            return {method_id, {option_id, false}};
        }
        break;
    case methods::MethodOptionType::integer: {
        int parsed = 0;
        const auto [end, error] =
            std::from_chars(value.data(), value.data() + value.size(), parsed);
        if (error == std::errc{} && end == value.data() + value.size()) {
            return {method_id, {option_id, parsed}};
        }
        break;
    }
    case methods::MethodOptionType::floating_point: {
        try {
            std::size_t end = 0;
            const auto parsed = std::stod(value, &end);
            if (end == value.size()) {
                return {method_id, {option_id, parsed}};
            }
        } catch (const std::exception&) {
            break;
        }
        break;
    }
    case methods::MethodOptionType::string:
        return {method_id, {option_id, value}};
    }
    throw std::runtime_error{"Invalid value for method option '" + method_id + "." + option_id +
                             "': " + value};
}

auto read_collection(const std::string& input_path,
                     const adapters::gemmi::InputOptions structural_options,
                     const bool structural_options_requested) -> ImportedCollection {
    std::ifstream input{input_path};
    if (!input) {
        throw std::runtime_error{"Unable to open input file: " + input_path};
    }

    auto extension = std::filesystem::path{input_path}.extension().string();
    std::ranges::transform(extension, extension.begin(), [](const unsigned char character) -> char {
        return static_cast<char>(std::tolower(character));
    });
    if (extension == ".sdf") {
        if (structural_options_requested) {
            throw std::runtime_error{
                "Structural input options are only supported for PDB and mmCIF input"};
        }
        adapters::native::sdf_input::SdfReader reader{input, input_path};
        return read_collection(reader, input_path, ImportedExportContext::Format::sdf);
    }
    if (extension == ".mol") {
        if (structural_options_requested) {
            throw std::runtime_error{
                "Structural input options are only supported for PDB and mmCIF input"};
        }
        adapters::native::mol_input::MolReader reader{input, input_path};
        return read_collection(reader, input_path, ImportedExportContext::Format::mol);
    }
    if (extension == ".mol2") {
        if (structural_options_requested) {
            throw std::runtime_error{
                "Structural input options are only supported for PDB and mmCIF input"};
        }
        adapters::native::mol2_input::Mol2Reader reader{input, input_path};
        return read_collection(reader, input_path, ImportedExportContext::Format::mol2);
    }
    if (extension == ".json") {
        if (structural_options_requested) {
            throw std::runtime_error{
                "Structural input options are only supported for PDB and mmCIF input"};
        }
        adapters::native::json_input::JsonReader reader{input, input_path,
                                                        structural_options.conformers};
        return read_collection(reader, input_path, ImportedExportContext::Format::json);
    }
    if (extension == ".pdb") {
        adapters::gemmi::pdb_input::PdbReader reader{input, input_path, structural_options};
        auto result = read_collection(reader, input_path, ImportedExportContext::Format::pdb);
        result.export_context.pdb_source = adapters::gemmi::mmcif_output::PdbSource{
            .structure = reader.source_structure(), .selection = reader.options().selection};
        return result;
    }
    if (extension == ".cif" || extension == ".mmcif") {
        adapters::gemmi::mmcif_input::MmcifReader reader{input, input_path, structural_options};
        auto result = read_collection(reader, input_path, ImportedExportContext::Format::mmcif);
        result.export_context.mmcif_source = adapters::gemmi::mmcif_output::MmcifSource{
            .document = reader.source_document(),
            .block_indices = reader.source_block_indices(),
            .selection = reader.options().selection};
        return result;
    }

    throw std::runtime_error{"Unsupported input file type: " + extension +
                             ". Supported types: .sdf, .mol, .mol2, .json, .pdb, .cif, .mmcif"};
}

} // namespace

void add_input_options(CLI::App& command, InputArguments& arguments) {
    command.add_option("input", arguments.path, "Input molecular file")->required();
    arguments.structural_selection_option =
        command.add_option("--structural-selection", arguments.structural_selection,
                           "PDB/mmCIF record selection: all, polymers-and-ligands, or polymers");
    arguments.structural_bonds_option =
        command.add_option("--structural-bonds", arguments.structural_bonds,
                           "PDB/mmCIF connectivity: none, explicit, templates, or hybrid");
    command.add_option("--conformers", arguments.conformer_selection,
                       "Conformers to read: first or all");
}

void add_selection_options(CLI::App& command, SelectionArguments& arguments) {
    arguments.method_option = command.add_option("--method", arguments.method_id, "Method ID");
    arguments.parameter_set_option =
        command.add_option("--parameter-set", arguments.parameter_set_id, "Parameter-set ID");
    command.add_flag("--permissive-types", arguments.permissive_types,
                     "Allow permissive parameter type classification");
    command.add_option("--method-option", arguments.method_options,
                       "Method option METHOD.OPTION=VALUE (repeatable)");
    command.add_option("--execution", arguments.execution,
                       "Execution: auto, full, cutoff, or cover");
    command.add_option("--radius", arguments.radius, "Cutoff or cover radius in angstrom");
    arguments.charge_correction_option =
        command.add_option("--charge-correction", arguments.charge_correction,
                           "Reduced-execution charge correction: uniform or none");
    arguments.cutoff_atom_threshold_option =
        command.add_option("--cutoff-atom-threshold", arguments.cutoff_atom_threshold,
                           "Automatic full-to-cutoff atom threshold, or unlimited");
    arguments.cover_atom_threshold_option =
        command.add_option("--cover-atom-threshold", arguments.cover_atom_threshold,
                           "Automatic cutoff-to-cover atom threshold, or unlimited");
    command.add_option("--threads", arguments.max_threads,
                       "Maximum calculation threads; 0 uses the oneTBB default");
}

auto import_input(const InputArguments& arguments) -> ImportedCollection {
    const auto options = adapters::gemmi::InputOptions{
        .selection = parse_record_selection(arguments.structural_selection),
        .bond_strategy = parse_bond_strategy(arguments.structural_bonds),
        .conformers = parse_conformer_selection(arguments.conformer_selection)};
    auto imported = read_collection(arguments.path, options,
                                    arguments.structural_selection_option->count() > 0 ||
                                        arguments.structural_bonds_option->count() > 0);
    imported.export_context.conformer_selection = arguments.conformer_selection;
    if (imported.export_context.format == ImportedExportContext::Format::pdb ||
        imported.export_context.format == ImportedExportContext::Format::mmcif) {
        imported.export_context.structural_input_policy =
            ImportedExportContext::StructuralInputPolicy{
                .selection = arguments.structural_selection, .bonds = arguments.structural_bonds};
    }
    return imported;
}

auto make_request(core::MoleculeCollection molecules, const SelectionArguments& arguments)
    -> calculation::AssessmentRequest {
    auto method_options = std::unordered_map<std::string, methods::MethodOptions>{};
    for (const auto& text : arguments.method_options) {
        const auto [method_id, option] = parse_method_option(text);
        method_options[method_id].set(option.first, option.second);
    }
    return {
        .molecules = std::move(molecules),
        .parameter_sets = parameters::load_default_parameter_sets(),
        .method_id = arguments.method_option->count() == 0 ? std::nullopt
                                                           : std::optional{arguments.method_id},
        .parameter_set_id = arguments.parameter_set_option->count() == 0
                                ? std::nullopt
                                : std::optional{arguments.parameter_set_id},
        .method_options = std::move(method_options),
        .classification_options = {.permissive_types = arguments.permissive_types},
        .execution_selection =
            calculation::ExecutionSelection{
                calculation::execution_selection_kind_from_string(arguments.execution),
                arguments.radius,
                arguments.charge_correction_option->count() == 0
                    ? std::nullopt
                    : std::optional{calculation::charge_correction_policy_from_string(
                          arguments.charge_correction)}},
        .resource_policy = {
            .cutoff_atom_threshold =
                arguments.cutoff_atom_threshold_option->count() == 0
                    ? std::optional<std::size_t>{calculation::default_cutoff_atom_threshold}
                    : parse_atom_threshold(arguments.cutoff_atom_threshold,
                                           "Cutoff atom threshold"),
            .cover_atom_threshold =
                arguments.cover_atom_threshold_option->count() == 0
                    ? std::optional<std::size_t>{calculation::default_cover_atom_threshold}
                    : parse_atom_threshold(arguments.cover_atom_threshold, "Cover atom threshold"),
            .max_threads = arguments.max_threads}};
}

} // namespace chargefw::cli
