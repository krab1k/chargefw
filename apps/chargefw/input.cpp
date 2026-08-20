#include "cli_support.h"

#include <chargefw/parameters/io/parameter_set_io.h>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <stdexcept>
#include <utility>

namespace chargefw::cli {
namespace {

template <typename Reader>
auto read_collection(Reader& reader, const std::string& input_path,
                     const ImportedCollection::Format format) -> ImportedCollection {
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
                              .records = std::move(records),
                              .format = format};
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

[[nodiscard]] auto parse_full_atom_threshold(const std::string& value)
    -> std::optional<std::size_t> {
    if (value == "unlimited") {
        return std::nullopt;
    }

    std::size_t threshold = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), threshold);
    if (error != std::errc{} || end != value.data() + value.size()) {
        throw std::runtime_error{
            "Full atom threshold must be a non-negative integer or 'unlimited'"};
    }
    return threshold;
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
        return read_collection(reader, input_path, ImportedCollection::Format::sdf);
    }
    if (extension == ".mol") {
        if (structural_options_requested) {
            throw std::runtime_error{
                "Structural input options are only supported for PDB and mmCIF input"};
        }
        adapters::native::mol_input::MolReader reader{input, input_path};
        return read_collection(reader, input_path, ImportedCollection::Format::mol);
    }
    if (extension == ".mol2") {
        if (structural_options_requested) {
            throw std::runtime_error{
                "Structural input options are only supported for PDB and mmCIF input"};
        }
        adapters::native::mol2_input::Mol2Reader reader{input, input_path};
        return read_collection(reader, input_path, ImportedCollection::Format::mol2);
    }
    if (extension == ".json") {
        if (structural_options_requested) {
            throw std::runtime_error{
                "Structural input options are only supported for PDB and mmCIF input"};
        }
        adapters::native::json_input::JsonReader reader{input, input_path,
                                                        structural_options.conformers};
        return read_collection(reader, input_path, ImportedCollection::Format::json);
    }
    if (extension == ".pdb") {
        adapters::gemmi::pdb_input::PdbReader reader{input, input_path, structural_options};
        auto result = read_collection(reader, input_path, ImportedCollection::Format::pdb);
        result.pdb_source = adapters::gemmi::mmcif_output::PdbSource{
            .structure = reader.source_structure(), .selection = reader.options().selection};
        return result;
    }
    if (extension == ".cif" || extension == ".mmcif") {
        adapters::gemmi::mmcif_input::MmcifReader reader{input, input_path, structural_options};
        auto result = read_collection(reader, input_path, ImportedCollection::Format::mmcif);
        result.mmcif_source = adapters::gemmi::mmcif_output::MmcifSource{
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
    command.add_option("--execution", arguments.execution,
                       "Execution: auto, full, cutoff, or cover");
    command.add_option("--radius", arguments.radius, "Cutoff or cover radius in angstrom");
    arguments.charge_correction_option =
        command.add_option("--charge-correction", arguments.charge_correction,
                           "Reduced-execution charge correction: uniform or none");
    arguments.full_atom_threshold_option =
        command.add_option("--full-atom-threshold", arguments.full_atom_threshold,
                           "Full-execution atom threshold, or unlimited");
}

auto import_input(const InputArguments& arguments) -> ImportedCollection {
    const auto options = adapters::gemmi::InputOptions{
        .selection = parse_record_selection(arguments.structural_selection),
        .bond_strategy = parse_bond_strategy(arguments.structural_bonds),
        .conformers = parse_conformer_selection(arguments.conformer_selection)};
    auto imported = read_collection(arguments.path, options,
                                    arguments.structural_selection_option->count() > 0 ||
                                        arguments.structural_bonds_option->count() > 0);
    imported.conformer_selection = arguments.conformer_selection;
    if (imported.format == ImportedCollection::Format::pdb ||
        imported.format == ImportedCollection::Format::mmcif) {
        imported.structural_input_policy = ImportedCollection::StructuralInputPolicy{
            .selection = arguments.structural_selection, .bonds = arguments.structural_bonds};
    }
    return imported;
}

auto make_request(const ImportedCollection& imported, const SelectionArguments& arguments)
    -> calculation::ApplicationCalculationRequest {
    return {.molecules = imported.molecules,
            .parameter_sets = parameters::load_default_parameter_sets(),
            .method_id = arguments.method_option->count() == 0 ? std::nullopt
                                                               : std::optional{arguments.method_id},
            .parameter_set_id = arguments.parameter_set_option->count() == 0
                                    ? std::nullopt
                                    : std::optional{arguments.parameter_set_id},
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
                .full_atom_threshold =
                    arguments.full_atom_threshold_option->count() == 0
                        ? std::optional<std::size_t>{calculation::default_full_atom_threshold}
                        : parse_full_atom_threshold(arguments.full_atom_threshold)}};
}

} // namespace chargefw::cli
