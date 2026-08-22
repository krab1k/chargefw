#include "cli_support.h"

#include <chargefw/charges/charge_collection.h>

#include <array>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <print>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <sys/resource.h>

namespace chargefw::cli {
namespace {

[[nodiscard]] auto result_document(const ImportedCollection& imported,
                                   const calculation::ApplicationAssessmentRequest& request,
                                   const calculation::ApplicationExecutionResult& result,
                                   const adapters::ExecutionMetrics& metrics)
    -> adapters::ChargeResultDocument;

[[nodiscard]] auto assignments_by_molecule(const charges::ChargeSet& charge_set,
                                           const std::size_t molecule_count)
    -> std::vector<charges::ChargeAssignment> {
    auto result = std::vector<charges::ChargeAssignment>{};
    result.reserve(molecule_count);
    for (std::size_t molecule_index = 0; molecule_index < molecule_count; ++molecule_index) {
        const auto found =
            std::ranges::find_if(charge_set.assignments(),
                                 [molecule_index](const charges::ChargeAssignment& assignment) {
                                     return assignment.target.molecule_index == molecule_index;
                                 });
        if (found == charge_set.assignments().end()) {
            throw std::runtime_error{"No charge assignment for molecule " +
                                     std::to_string(molecule_index)};
        }
        const auto remaining = charge_set.assignments().subspan(
            static_cast<std::size_t>(std::distance(charge_set.assignments().begin(), found)) + 1);
        const auto duplicate = std::ranges::find_if(
            remaining, [molecule_index](const charges::ChargeAssignment& assignment) {
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

void write_json(const std::filesystem::path& path, const adapters::ChargeResultDocument& document) {
    auto output = std::ofstream{path};
    if (!output) {
        throw std::runtime_error{"Unable to open output file: " + path.string()};
    }
    adapters::native::json_output::JsonWriter{output}.write(document);
}

void write_mmcif(const std::filesystem::path& path, const ImportedCollection& imported,
                 const charges::ChargeSet& charges) {
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

void write_mol2(const std::filesystem::path& path, const std::string& input_path,
                const ImportedCollection& imported,
                const std::span<const charges::ChargeAssignment> assignments) {
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

void write_sdf(const std::filesystem::path& path, const std::string& input_path,
               const ImportedCollection& imported,
               const std::span<const charges::ChargeAssignment> assignments,
               const std::string_view method) {
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
                                   const calculation::ApplicationAssessmentRequest& request,
                                   const calculation::ApplicationExecutionResult& result,
                                   const adapters::ExecutionMetrics& metrics)
    -> adapters::ChargeResultDocument {
    auto warnings = std::vector<std::string>{};
    warnings.reserve(result.execution_issues.size());
    for (const auto& issue : result.execution_issues) {
        warnings.push_back(issue.message);
    }
    auto provenance = adapters::CalculationProvenance{
        .requested = {.method_id = request.method_id,
                      .parameter_set_id = request.parameter_set_id,
                      .permissive_types = request.classification_options.permissive_types,
                      .full_atom_threshold = request.resource_policy.full_atom_threshold,
                      .max_threads = request.resource_policy.max_threads,
                      .execution_kind =
                          std::string{calculation::to_string(request.execution_selection.kind())},
                      .execution_radius = request.execution_selection.radius(),
                      .execution_charge_correction =
                          request.execution_selection.charge_correction().transform(
                              [](const calculation::ChargeCorrectionPolicy policy) {
                                  return std::string{calculation::to_string(policy)};
                              }),
                      .structural_input_policy =
                          imported.structural_input_policy.has_value()
                              ? std::optional{adapters::StructuralInputPolicyProvenance{
                                    .selection = imported.structural_input_policy->selection,
                                    .bonds = imported.structural_input_policy->bonds}}
                              : std::nullopt,
                      .conformer_selection = imported.conformer_selection},
        .effective = {.method_id = result.charges.has_value()
                                       ? std::optional{std::string{result.charges->method_id()}}
                                       : std::nullopt,
                      .parameter_set_id =
                          result.charges.has_value()
                              ? result.charges->parameter_set_id().transform(
                                    [](const std::string_view id) { return std::string{id}; })
                              : std::nullopt,
                      .execution_mode = result.execution_policy.has_value()
                                            ? std::optional{std::string{calculation::to_string(
                                                  result.execution_policy->mode())}}
                                            : std::nullopt,
                      .execution_radius = result.execution_policy.has_value()
                                              ? result.execution_policy->radius()
                                              : std::nullopt,
                      .execution_charge_correction =
                          result.execution_policy.has_value()
                              ? std::optional{std::string{calculation::to_string(
                                    result.execution_policy->charge_correction())}}
                              : std::nullopt,
                      .warnings = std::move(warnings)},
        .execution_metrics = metrics};
    for (const auto& [method_id, options] : request.method_options) {
        provenance.requested.method_options.emplace(method_id, options);
    }
    if (result.effective_method_options.has_value() && result.charges.has_value()) {
        provenance.effective.method_options.emplace(std::string{result.charges->method_id()},
                                                    *result.effective_method_options);
    }
    auto document = adapters::ChargeResultDocument{.generator_name = "ChargeFW",
                                                   .generator_version = "0.0.1",
                                                   .records = {},
                                                   .calculation_provenance = provenance};
    document.records.reserve(imported.molecules.size());
    for (std::size_t index = 0; index < imported.molecules.size(); ++index) {
        document.records.push_back(adapters::ChargeResultRecord{
            .identity = imported.records[index].identity, .charges = result.charges});
    }
    return document;
}

} // namespace

auto utc_timestamp() -> std::string {
    const auto now = std::chrono::system_clock::now();
    const auto time = std::chrono::system_clock::to_time_t(now);
    std::tm utc{};
    gmtime_r(&time, &utc);
    const auto milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
    auto output = std::ostringstream{};
    output << std::put_time(&utc, "%Y-%m-%dT%H:%M:%S") << '.' << std::setfill('0') << std::setw(3)
           << milliseconds.count() << 'Z';
    return output.str();
}

auto peak_resident_memory_mb() -> double {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return 0.0;
    }
    return static_cast<double>(usage.ru_maxrss) / 1024.0;
}

auto write_calculation_outputs(const std::string& output_directory, const std::string& input_path,
                               const ImportedCollection& imported,
                               const calculation::ApplicationAssessmentRequest& request,
                               const calculation::ApplicationExecutionResult& result,
                               CalculationRun& run) -> int {
    run.metrics.peak_resident_memory_mb = peak_resident_memory_mb();
    if (!result.calculated()) {
        run.metrics.ended_at = utc_timestamp();
        run.metrics.runtime_seconds =
            std::chrono::duration<double>{std::chrono::steady_clock::now() - run.started}.count();
        adapters::native::json_output::JsonWriter{std::cout}.write(
            result_document(imported, request, result, run.metrics));
        return 1;
    }

    const auto structural_output = imported.format == ImportedCollection::Format::pdb ||
                                   imported.format == ImportedCollection::Format::mmcif;
    if (!structural_output && imported.format == ImportedCollection::Format::json &&
        std::ranges::any_of(imported.molecules.molecules(), [](const core::Molecule& molecule) {
            return molecule.conformer_count() > 1;
        })) {
        throw std::runtime_error{
            "JSON input with multiple conformers cannot be written to SDF or MOL2"};
    }

    const auto directory = std::filesystem::path{output_directory};
    std::error_code directory_error;
    std::filesystem::create_directories(directory, directory_error);
    if (directory_error) {
        throw std::runtime_error{"Unable to create output directory: " + directory.string() + ": " +
                                 directory_error.message()};
    }
    if (!std::filesystem::is_directory(directory)) {
        throw std::runtime_error{"Output path is not a directory: " + directory.string()};
    }
    const auto prefix =
        directory / (std::filesystem::path{input_path}.stem().string() + ".chargefw");
    const auto* charges = result.charges ? std::addressof(*result.charges) : nullptr;
    if (charges == nullptr) {
        throw std::runtime_error{"calculation result is missing charges"};
    }
    const auto writing_started = std::chrono::steady_clock::now();
    write_mmcif(prefix.string() + ".cif", imported, *charges);
    if (structural_output) {
        std::println("Wrote {} and {}", prefix.string() + ".json", prefix.string() + ".cif");
    } else {
        const auto assignments = assignments_by_molecule(*charges, imported.molecules.size());
        write_sdf(prefix.string() + ".sdf", input_path, imported, assignments,
                  charges->method_id());
        write_mol2(prefix.string() + ".mol2", input_path, imported, assignments);
        std::println("Wrote {}, {}, {}, and {}", prefix.string() + ".json",
                     prefix.string() + ".sdf", prefix.string() + ".mol2", prefix.string() + ".cif");
    }
    run.metrics.writing_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - writing_started}.count();
    run.metrics.peak_resident_memory_mb = peak_resident_memory_mb();
    run.metrics.ended_at = utc_timestamp();
    run.metrics.runtime_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - run.started}.count();
    write_json(prefix.string() + ".json", result_document(imported, request, result, run.metrics));
    return 0;
}

} // namespace chargefw::cli
