#include <CLI/CLI.hpp>
#include <algorithm>
#include <chargefw/adapters/gemmi/all.h>
#include <chargefw/adapters/native/all.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <array>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
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
        CLI::App app{"Calculate empirical partial atomic charges from a molecular file."};
        std::string input_path;
        std::string output_directory;
        std::string structural_selection = "all";
        std::string structural_bonds = "hybrid";
        app.add_option("input", input_path,
                       "Input .sdf, .mol, .mol2, .pdb, .cif, .mmcif, or ChargeFW .json file")
            ->required();
        app.add_option("output", output_directory, "Output directory")->required();
        const auto* structural_selection_option =
            app.add_option("--structural-selection", structural_selection,
                           "PDB/mmCIF record selection: all, polymers-and-ligands, or polymers")
                ->default_val("all");
        const auto* structural_bonds_option =
            app.add_option("--structural-bonds", structural_bonds,
                           "PDB/mmCIF connectivity: none, explicit, templates, or hybrid")
                ->default_val("hybrid");
        CLI11_PARSE(app, argc, argv);

        const auto structural_options =
            adapters::gemmi::InputOptions{.selection = parse_record_selection(structural_selection),
                                          .bond_strategy = parse_bond_strategy(structural_bonds)};
        const auto structural_options_requested =
            structural_selection_option->count() > 0 || structural_bonds_option->count() > 0;
        const auto imported =
            read_collection(input_path, structural_options, structural_options_requested);
        const auto parameter_sets = parameters::load_default_parameter_sets();

        const auto result = calculation::calculate(calculation::ApplicationCalculationRequest{
            .molecules = imported.molecules, .parameter_sets = parameter_sets});

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
        const auto input_name = std::filesystem::path{input_path}.stem();
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
        write_sdf(prefix.string() + ".sdf", input_path, imported, assignments,
                  charges->method_id());
        write_mol2(prefix.string() + ".mol2", input_path, imported, assignments);
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
