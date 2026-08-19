#include <CLI/CLI.hpp>
#include <algorithm>
#include <chargefw/adapters/gemmi/all.h>
#include <chargefw/adapters/native/all.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <array>
#include <cctype>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <print>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;
namespace adapters = chargefw::adapters;

namespace {

struct ImportedCollection {
    enum class Format : std::uint8_t { sdf, mol, mol2, json, pdb, mmcif };

    core::MoleculeCollection molecules;
    std::vector<adapters::ImportedMoleculeRecord> records;
    std::optional<adapters::gemmi::mmcif_output::MmcifSource> mmcif_source;
    std::optional<adapters::gemmi::mmcif_output::PdbSource> pdb_source;
    Format format;
};

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

[[nodiscard]] auto parse_execution_selection(const std::string& value)
    -> calculation::ExecutionSelectionKind {
    if (value == "auto") {
        return calculation::ExecutionSelectionKind::automatic;
    }
    if (value == "full") {
        return calculation::ExecutionSelectionKind::full;
    }
    if (value == "cutoff") {
        return calculation::ExecutionSelectionKind::cutoff;
    }
    if (value == "cover") {
        return calculation::ExecutionSelectionKind::cover;
    }
    throw std::runtime_error{"Unsupported execution selection: " + value};
}

[[nodiscard]] auto parse_charge_correction(const std::string& value)
    -> calculation::ChargeCorrectionPolicy {
    if (value == "none") {
        return calculation::ChargeCorrectionPolicy::none;
    }
    if (value == "uniform") {
        return calculation::ChargeCorrectionPolicy::uniform;
    }
    throw std::runtime_error{"Unsupported charge correction policy: " + value};
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
        adapters::native::json_input::JsonReader reader{input, input_path};
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

[[nodiscard]] auto assignments_by_molecule(const charges::ChargeSet& charge_set,
                                           const std::size_t molecule_count)
    -> std::vector<charges::ChargeAssignment> {
    auto result = std::vector<charges::ChargeAssignment>{};
    result.reserve(molecule_count);
    for (std::size_t molecule_index = 0; molecule_index < molecule_count; ++molecule_index) {
        const auto found = std::ranges::find_if(
            charge_set.assignments(),
            [molecule_index](const charges::ChargeAssignment& assignment) -> bool {
                return assignment.target.molecule_index == molecule_index;
            });
        if (found == charge_set.assignments().end()) {
            throw std::runtime_error{"No charge assignment for molecule " +
                                     std::to_string(molecule_index)};
        }
        const auto remaining = charge_set.assignments().subspan(
            static_cast<std::size_t>(std::distance(charge_set.assignments().begin(), found)) + 1);
        const auto duplicate = std::ranges::find_if(
            remaining, [molecule_index](const charges::ChargeAssignment& assignment) -> bool {
                return assignment.target.molecule_index == molecule_index;
            });
        if (duplicate != remaining.end()) {
            throw std::runtime_error{
                "Molecular output does not support multiple conformer assignments per molecule"};
        }
        result.push_back(*found);
    }
    return result;
}

auto write_json(const std::filesystem::path& path, const adapters::ChargeResultDocument& document)
    -> void {
    auto output = std::ofstream{path};
    if (!output) {
        throw std::runtime_error{"Unable to open output file: " + path.string()};
    }
    adapters::native::json_output::JsonWriter{output}.write(document);
}

auto write_mmcif(const std::filesystem::path& path, const ImportedCollection& imported,
                 const charges::ChargeSet& charges) -> void {
    auto output = std::ofstream{path};
    if (!output) {
        throw std::runtime_error{"Unable to open output file: " + path.string()};
    }
    auto writer = adapters::gemmi::mmcif_output::MmcifWriter{output};
    if (imported.mmcif_source.has_value()) {
        writer.write_mmcif(imported.records, charges, *imported.mmcif_source, "ChargeFW", "0.0.1");
    } else if (imported.pdb_source.has_value()) {
        writer.write_pdb(imported.records.front(), charges, *imported.pdb_source, "ChargeFW",
                         "0.0.1");
    } else {
        writer.write_generated(imported.records, charges, "ChargeFW", "0.0.1");
    }
}

auto write_mol2(const std::filesystem::path& path, const std::string& input_path,
                const ImportedCollection& imported,
                const std::span<const charges::ChargeAssignment> assignments) -> void {
    auto output = std::ofstream{path, std::ios::binary};
    if (!output) {
        throw std::runtime_error{"Unable to open output file: " + path.string()};
    }
    auto writer = adapters::native::mol2_output::Mol2Writer{output};
    if (imported.format == ImportedCollection::Format::mol2) {
        writer.write_preserving_source(input_path, assignments);
        return;
    }
    for (std::size_t index = 0; index < imported.molecules.size(); ++index) {
        writer.write_generated(imported.molecules[index], assignments[index]);
    }
}

auto write_sdf(const std::filesystem::path& path, const std::string& input_path,
               const ImportedCollection& imported,
               const std::span<const charges::ChargeAssignment> assignments,
               const std::string_view method) -> void {
    auto output = std::ofstream{path, std::ios::binary};
    if (!output) {
        throw std::runtime_error{"Unable to open output file: " + path.string()};
    }
    auto writer = adapters::native::sdf_output::SdfWriter{output};
    const auto properties = std::array{adapters::native::sdf_output::ChargeProperty{
        .charge_type_id = 1, .assignments = assignments, .method = method}};
    if (imported.format == ImportedCollection::Format::sdf) {
        writer.write_preserving_source(input_path, properties);
        return;
    }
    for (std::size_t index = 0; index < imported.molecules.size(); ++index) {
        const auto property = std::array{adapters::native::sdf_output::ChargeProperty{
            .charge_type_id = 1, .assignments = assignments.subspan(index, 1), .method = method}};
        writer.write_generated(imported.molecules[index], property,
                               adapters::native::sdf_output::MolFormat::v2000);
    }
}

struct InputArguments {
    std::string path;
    std::string structural_selection = "all";
    std::string structural_bonds = "hybrid";
    CLI::Option* structural_selection_option = nullptr;
    CLI::Option* structural_bonds_option = nullptr;
};

struct SelectionArguments {
    std::string method_id;
    std::string parameter_set_id;
    std::string execution = "auto";
    std::optional<double> radius;
    std::string charge_correction;
    std::string full_atom_threshold;
    bool permissive_types = false;
    CLI::Option* method_option = nullptr;
    CLI::Option* parameter_set_option = nullptr;
    CLI::Option* charge_correction_option = nullptr;
    CLI::Option* full_atom_threshold_option = nullptr;
};

auto add_input_options(CLI::App& command, InputArguments& arguments) -> void {
    command.add_option("input", arguments.path, "Input molecular file")->required();
    arguments.structural_selection_option =
        command.add_option("--structural-selection", arguments.structural_selection,
                           "PDB/mmCIF record selection: all, polymers-and-ligands, or polymers");
    arguments.structural_bonds_option =
        command.add_option("--structural-bonds", arguments.structural_bonds,
                           "PDB/mmCIF connectivity: none, explicit, templates, or hybrid");
}

auto add_selection_options(CLI::App& command, SelectionArguments& arguments) -> void {
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

[[nodiscard]] auto import_input(const InputArguments& arguments) -> ImportedCollection {
    const auto options = adapters::gemmi::InputOptions{
        .selection = parse_record_selection(arguments.structural_selection),
        .bond_strategy = parse_bond_strategy(arguments.structural_bonds)};
    return read_collection(arguments.path, options,
                           arguments.structural_selection_option->count() > 0 ||
                               arguments.structural_bonds_option->count() > 0);
}

[[nodiscard]] auto make_request(const ImportedCollection& imported,
                                const SelectionArguments& arguments)
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
                    parse_execution_selection(arguments.execution), arguments.radius,
                    arguments.charge_correction_option->count() == 0
                        ? std::nullopt
                        : std::optional{parse_charge_correction(arguments.charge_correction)}},
            .resource_policy = {
                .full_atom_threshold =
                    arguments.full_atom_threshold_option->count() == 0
                        ? std::optional<std::size_t>{calculation::default_full_atom_threshold}
                        : parse_full_atom_threshold(arguments.full_atom_threshold)}};
}

auto print_inspection(const ImportedCollection& imported) -> void {
    std::println("records: {}", imported.molecules.size());
    for (std::size_t index = 0; index < imported.molecules.size(); ++index) {
        const auto& molecule = imported.molecules[index];
        std::map<int, std::size_t> elements;
        for (const auto& atom : molecule.atoms()) {
            ++elements[atom.atomic_number()];
        }
        std::print("record {} ({}) atoms={} bonds={} conformers={} coordinates={} formal_charge={} "
                   "elements=",
                   index, imported.records[index].identity.record_id, molecule.atom_count(),
                   molecule.bond_count(), molecule.conformer_count(), molecule.has_coordinates(),
                   core::total_formal_charge(molecule));
        bool first = true;
        for (const auto& [atomic_number, count] : elements) {
            std::print("{}{}:{}", first ? "" : ",", core::element_symbol(atomic_number), count);
            first = false;
        }
        std::println();
    }
}

auto print_applicability(const calculation::ApplicationAssessmentResult& assessment) -> void {
    std::println("applicable candidates: {}", assessment.applicability.applicable.size());
    for (const auto& candidate : assessment.applicability.applicable) {
        std::print("applicable method={} parameter_set={}", candidate.method->id(),
                   candidate.parameter_set == nullptr ? "-" : candidate.parameter_set->id());
        for (const auto& execution : candidate.execution_assessments) {
            const auto mode = execution.mode == calculation::ExecutionMode::full     ? "full"
                              : execution.mode == calculation::ExecutionMode::cutoff ? "cutoff"
                                                                                     : "cover";
            std::print(" {}={}", mode,
                       execution.availability == methods::ExecutionAvailability::unsupported
                           ? "unsupported"
                           : "available");
        }
        std::println();
    }
    std::println("rejected candidates: {}", assessment.applicability.rejected.size());
    const auto& registry = methods::method_registry();
    for (const auto& rejected : assessment.applicability.rejected) {
        std::print("rejected method={}", registry.methods()[rejected.method_index]->id());
        if (rejected.parameter_set_index.has_value()) {
            std::print(" parameter_set={}",
                       assessment.parameter_sets[*rejected.parameter_set_index].id());
        }
        for (const auto& issue : rejected.issues) {
            std::print("; {}", issue.message);
        }
        std::println();
    }
    if (assessment.selected == nullptr) {
        std::println("selected execution: none");
        return;
    }
    std::println(
        "selected method={} parameter_set={} execution={}", assessment.selected->method->id(),
        assessment.selected->parameter_set == nullptr ? "-"
                                                      : assessment.selected->parameter_set->id(),
        assessment.execution_policy->mode() == calculation::ExecutionMode::full     ? "full"
        : assessment.execution_policy->mode() == calculation::ExecutionMode::cutoff ? "cutoff"
                                                                                    : "cover");
}

[[nodiscard]] auto result_document(const ImportedCollection& imported,
                                   const calculation::ApplicationCalculationResult& result)
    -> adapters::ChargeResultDocument {
    auto document = adapters::ChargeResultDocument{
        .generator_name = "ChargeFW", .generator_version = "0.0.1", .records = {}};
    document.records.reserve(imported.molecules.size());
    for (std::size_t index = 0; index < imported.molecules.size(); ++index) {
        document.records.push_back(
            adapters::ChargeResultRecord{.identity = imported.records[index].identity,
                                         .mapping = imported.records[index].mapping,
                                         .charges = result.charges});
    }
    return document;
}

} // namespace

auto run(int argc, char* argv[]) -> int {
    try {
        CLI::App app{"ChargeFW molecular charge calculation and inspection."};
        InputArguments calculate_input;
        InputArguments inspect_input;
        InputArguments applicability_input;
        SelectionArguments calculate_selection;
        SelectionArguments applicability_selection;
        std::string output_directory;
        std::string parameter_method;
        auto* calculate = app.add_subcommand("calculate", "Calculate and write partial charges");
        auto* inspect = app.add_subcommand("inspect", "Inspect imported molecular records");
        auto* applicability =
            app.add_subcommand("applicability", "Report applicable charge methods");
        auto* methods_command = app.add_subcommand("methods", "List registered charge methods");
        auto* parameters_command = app.add_subcommand("parameters", "List bundled parameter sets");
        add_input_options(*calculate, calculate_input);
        calculate->add_option("output", output_directory, "Output directory")->required();
        add_selection_options(*calculate, calculate_selection);
        add_input_options(*inspect, inspect_input);
        add_input_options(*applicability, applicability_input);
        add_selection_options(*applicability, applicability_selection);
        parameters_command->add_option("method", parameter_method, "Limit results to a method ID");
        CLI11_PARSE(app, argc, argv);

        if (*methods_command) {
            for (const auto& method : methods::method_registry().methods()) {
                std::println("{}\t{}", method->id(), method->metadata().name);
            }
            return 0;
        }
        if (*parameters_command) {
            for (const auto& parameter_set : parameters::load_default_parameter_sets()) {
                if (parameter_method.empty() || parameter_set.method_id() == parameter_method) {
                    std::println("{}\t{}\t{}", parameter_set.id(), parameter_set.method_id(),
                                 parameter_set.name());
                }
            }
            return 0;
        }
        if (*inspect) {
            print_inspection(import_input(inspect_input));
            return 0;
        }
        if (*applicability) {
            const auto imported = import_input(applicability_input);
            const auto request = make_request(imported, applicability_selection);
            print_applicability(calculation::assess(request));
            return 0;
        }
        if (!*calculate) {
            throw std::invalid_argument{"a subcommand is required; use calculate, inspect, "
                                        "applicability, methods, or parameters"};
        }

        const auto imported = import_input(calculate_input);
        const auto request = make_request(imported, calculate_selection);
        const auto result = calculation::calculate(request);

        if (!result.calculated()) {
            adapters::native::json_output::JsonWriter{std::cout}.write(
                result_document(imported, result));
            return 1;
        }

        const auto structural_output = imported.format == ImportedCollection::Format::pdb ||
                                       imported.format == ImportedCollection::Format::mmcif;
        if (!structural_output && imported.format == ImportedCollection::Format::json &&
            std::ranges::any_of(imported.molecules.molecules(),
                                [](const core::Molecule& molecule) -> bool {
                                    return molecule.conformer_count() > 1;
                                })) {
            throw std::runtime_error{
                "JSON input with multiple conformers cannot be written to SDF or MOL2"};
        }

        const auto directory = std::filesystem::path{output_directory};
        std::error_code directory_error;
        std::filesystem::create_directories(directory, directory_error);
        if (directory_error) {
            throw std::runtime_error{"Unable to create output directory: " + directory.string() +
                                     ": " + directory_error.message()};
        }
        if (!std::filesystem::is_directory(directory)) {
            throw std::runtime_error{"Output path is not a directory: " + directory.string()};
        }
        const auto input_name = std::filesystem::path{calculate_input.path}.stem();
        const auto prefix = directory / (input_name.string() + ".chargefw");
        const auto output = result_document(imported, result);
        write_json(prefix.string() + ".json", output);
        const auto* charges = result.charges ? std::addressof(*result.charges) : nullptr;
        if (charges == nullptr) {
            throw std::runtime_error{"calculation result is missing charges"};
        }
        write_mmcif(prefix.string() + ".cif", imported, *charges);
        if (structural_output) {
            std::println("Wrote {} and {}", prefix.string() + ".json", prefix.string() + ".cif");
            return 0;
        }
        const auto assignments = assignments_by_molecule(*charges, imported.molecules.size());
        write_sdf(prefix.string() + ".sdf", calculate_input.path, imported, assignments,
                  charges->method_id());
        write_mol2(prefix.string() + ".mol2", calculate_input.path, imported, assignments);
        std::println("Wrote {}, {}, {}, and {}", prefix.string() + ".json",
                     prefix.string() + ".sdf", prefix.string() + ".mol2", prefix.string() + ".cif");

        return 0;
    } catch (const std::exception& error) {
        std::print(std::cerr, "Fatal error: {}\n", error.what());
        return 1;
    }
}

auto main(int argc, char* argv[]) noexcept -> int {
    try {
        return run(argc, argv);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "Fatal error: %s\n", error.what());
        return 1;
    }
}
