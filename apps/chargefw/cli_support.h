#pragma once

#include <CLI/CLI.hpp>
#include <chargefw/adapters/gemmi/all.h>
#include <chargefw/adapters/native/all.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/core/molecule_collection.h>

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace chargefw::cli {

struct ImportedCollection {
    enum class Format : std::uint8_t { sdf, mol, mol2, json, pdb, mmcif };

    struct StructuralInputPolicy {
        std::string selection;
        std::string bonds;
    };

    core::MoleculeCollection molecules;
    std::vector<adapters::ImportedMoleculeRecord> records;
    std::optional<adapters::gemmi::mmcif_output::MmcifSource> mmcif_source;
    std::optional<adapters::gemmi::mmcif_output::PdbSource> pdb_source;
    std::optional<StructuralInputPolicy> structural_input_policy;
    Format format;
};

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

void add_input_options(CLI::App& command, InputArguments& arguments);
void add_selection_options(CLI::App& command, SelectionArguments& arguments);
[[nodiscard]] auto import_input(const InputArguments& arguments) -> ImportedCollection;
[[nodiscard]] auto make_request(const ImportedCollection& imported,
                                const SelectionArguments& arguments)
    -> calculation::ApplicationCalculationRequest;

void print_inspection(const ImportedCollection& imported);
void print_applicability(const calculation::ApplicationAssessmentResult& assessment);
void print_methods();
void print_parameter_sets(const std::string& method_id);

[[nodiscard]] auto
write_calculation_outputs(const std::string& output_directory, const std::string& input_path,
                          const ImportedCollection& imported,
                          const calculation::ApplicationCalculationRequest& request,
                          const calculation::ApplicationCalculationResult& result) -> int;

} // namespace chargefw::cli
