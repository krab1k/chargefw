#include "cli_support.h"

#include <chargefw/charges/charge_collection.h>
#include <chargefw/config.h>

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
#include <set>
#include <sstream>
#include <stdexcept>
#include <sys/resource.h>
#include <variant>

namespace chargefw::cli {
namespace {

[[nodiscard]] auto result_document(const ImportedExportContext& export_context,
                                   const adapters::RequestedCalculationProvenance& requested,
                                   const calculation::ExecutionResult& result,
                                   const adapters::ExecutionMetrics& metrics)
    -> adapters::ChargeResultDocument;

[[nodiscard]] auto prerequisite_code(const methods::PrerequisiteIssueKind kind)
    -> std::string_view {
    switch (kind) {
    case methods::PrerequisiteIssueKind::invalid_options:
        return "invalid_method_options";
    case methods::PrerequisiteIssueKind::missing_feature:
        return "missing_feature";
    case methods::PrerequisiteIssueKind::invalid_geometry:
        return "invalid_geometry";
    case methods::PrerequisiteIssueKind::unsupported_molecule:
        return "unsupported_molecule";
    case methods::PrerequisiteIssueKind::missing_parameters:
        return "missing_parameters";
    case methods::PrerequisiteIssueKind::parameter_classification_failed:
        return "parameter_classification_failed";
    }
    throw std::logic_error{"unknown prerequisite issue kind"};
}

void append_unique(std::vector<adapters::ResultDiagnostic>& diagnostics,
                   adapters::ResultDiagnostic diagnostic) {
    const auto duplicate =
        std::ranges::any_of(diagnostics, [&diagnostic](const adapters::ResultDiagnostic& existing) {
            return existing.code == diagnostic.code && existing.message == diagnostic.message;
        });
    if (!duplicate) {
        diagnostics.push_back(std::move(diagnostic));
    }
}

void report_diagnostics(const adapters::ChargeResultDocument& document) {
    auto reported = std::set<std::pair<std::string, std::string>>{};
    const auto report = [&reported](const adapters::ResultDiagnostic& diagnostic) {
        if (!reported.emplace(diagnostic.code, diagnostic.message).second) {
            return;
        }
        const auto label = diagnostic.severity == adapters::DiagnosticSeverity::warning ? "Warning"
                           : diagnostic.severity == adapters::DiagnosticSeverity::info  ? "Info"
                                                                                        : "Error";
        std::println(std::cerr, "{}: {}", label, diagnostic.message);
    };
    for (const auto& diagnostic : document.diagnostics) {
        report(diagnostic);
    }
    for (const auto& record : document.records) {
        for (const auto& diagnostic : record.diagnostics) {
            report(diagnostic);
        }
    }
}

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

void write_mmcif(const std::filesystem::path& path, const ImportedExportContext& export_context,
                 const charges::ChargeSet& charges) {
    auto output = std::ofstream{path};
    if (!output) {
        throw std::runtime_error{"Unable to open output file: " + path.string()};
    }
    auto writer = adapters::gemmi::mmcif_output::MmcifWriter{output};
    if (export_context.mmcif_source.has_value()) {
        writer.write_mmcif(export_context.records, charges, *export_context.mmcif_source,
                           "ChargeFW", CHARGEFW_VERSION_STRING);
    } else if (export_context.pdb_source.has_value()) {
        writer.write_pdb(export_context.records.front(), charges, *export_context.pdb_source,
                         "ChargeFW", CHARGEFW_VERSION_STRING);
    } else {
        writer.write_generated(export_context.records, charges, "ChargeFW",
                               CHARGEFW_VERSION_STRING);
    }
}

void write_mol2(const std::filesystem::path& path, const std::string& input_path,
                const ImportedExportContext& export_context,
                const std::span<const charges::ChargeAssignment> assignments) {
    auto output = std::ofstream{path, std::ios::binary};
    if (!output) {
        throw std::runtime_error{"Unable to open output file: " + path.string()};
    }
    auto writer = adapters::native::mol2_output::Mol2Writer{output};
    if (export_context.format == ImportedExportContext::Format::mol2) {
        writer.write_preserving_source(input_path, assignments);
        return;
    }
    for (std::size_t index = 0; index < export_context.records.size(); ++index) {
        writer.write_generated(export_context.records[index].molecule, assignments[index]);
    }
}

void write_sdf(const std::filesystem::path& path, const std::string& input_path,
               const ImportedExportContext& export_context, const charges::ChargeSet& charge_set) {
    auto output = std::ofstream{path, std::ios::binary};
    if (!output) {
        throw std::runtime_error{"Unable to open output file: " + path.string()};
    }
    auto writer = adapters::native::sdf_output::SdfWriter{output};
    const auto assignments = charge_set.assignments();
    const auto properties = std::array{adapters::native::sdf_output::ChargeProperty{
        .charge_type_id = 1,
        .assignments = assignments,
        .method = charge_set.method_id(),
        .parameter_set = charge_set.parameter_set_id().value_or(""),
        .software_name = "ChargeFW",
        .software_version = CHARGEFW_VERSION_STRING}};
    if (export_context.format == ImportedExportContext::Format::sdf) {
        writer.write_preserving_source(input_path, properties);
        return;
    }
    for (std::size_t index = 0; index < export_context.records.size(); ++index) {
        const auto property = std::array{adapters::native::sdf_output::ChargeProperty{
            .charge_type_id = 1,
            .assignments = assignments.subspan(index, 1),
            .method = charge_set.method_id(),
            .parameter_set = charge_set.parameter_set_id().value_or(""),
            .software_name = "ChargeFW",
            .software_version = CHARGEFW_VERSION_STRING}};
        writer.write_generated(export_context.records[index].molecule, property,
                               adapters::native::sdf_output::MolFormat::v2000);
    }
}

[[nodiscard]] auto result_document(const ImportedExportContext& export_context,
                                   const adapters::RequestedCalculationProvenance& requested,
                                   const calculation::ExecutionResult& result,
                                   const adapters::ExecutionMetrics& metrics)
    -> adapters::ChargeResultDocument {
    const auto status = result.status;
    const auto diagnostic = [&result] -> std::optional<adapters::ResultDiagnostic> {
        switch (result.status) {
        case calculation::ExecutionStatus::success:
            return std::nullopt;
        case calculation::ExecutionStatus::invalid_input_or_request:
            return adapters::ResultDiagnostic{.severity = adapters::DiagnosticSeverity::error,
                                              .code = "invalid_input_or_request",
                                              .message = result.failure_message.value_or(
                                                  "Invalid input or calculation request.")};
        case calculation::ExecutionStatus::no_executable_plan:
            return adapters::ResultDiagnostic{
                .severity = adapters::DiagnosticSeverity::error,
                .code = "no_executable_plan",
                .message = "No executable method and parameter-set plan was found."};
        case calculation::ExecutionStatus::numerical_failure:
            return adapters::ResultDiagnostic{
                .severity = adapters::DiagnosticSeverity::error,
                .code = "numerical_failure",
                .message = result.failure_message.value_or("Calculation failed numerically.")};
        case calculation::ExecutionStatus::cancelled:
            return adapters::ResultDiagnostic{.severity = adapters::DiagnosticSeverity::info,
                                              .code = "calculation_cancelled",
                                              .message =
                                                  "Calculation was cancelled before completion."};
        }
        throw std::logic_error{"unknown execution result status"};
    }();
    auto effective = adapters::EffectiveCalculationProvenance{};
    if (result.effective.has_value()) {
        const auto& calculation = *result.effective;
        effective.method_id = calculation.method_id;
        effective.parameter_set_id = calculation.parameter_set_id;
        effective.execution_mode =
            std::string{calculation::to_string(calculation.execution_policy.mode())};
        effective.execution_radius = calculation.execution_policy.radius();
        effective.execution_charge_correction =
            std::string{calculation::to_string(calculation.execution_policy.charge_correction())};
        effective.warnings.reserve(calculation.execution_issues.size());
        for (const auto& issue : calculation.execution_issues) {
            effective.warnings.push_back(issue.message);
        }
        effective.method_options.emplace(calculation.method_id, calculation.method_options);
    }
    auto provenance = adapters::CalculationProvenance{
        .requested = requested, .effective = std::move(effective), .execution_metrics = metrics};
    auto document = adapters::ChargeResultDocument{
        .generator_name = "ChargeFW",
        .generator_version = CHARGEFW_VERSION_STRING,
        .status = status,
        .diagnostics = diagnostic.has_value() ? std::vector{*diagnostic}
                                              : std::vector<adapters::ResultDiagnostic>{},
        .records = {},
        .calculation_provenance = provenance};
    document.records.reserve(export_context.records.size());
    for (std::size_t molecule_index = 0; molecule_index < export_context.records.size();
         ++molecule_index) {
        const auto& record = export_context.records[molecule_index];
        auto record_diagnostics = std::vector<adapters::ResultDiagnostic>{};
        for (const auto& import_diagnostic : record.diagnostics) {
            append_unique(record_diagnostics, adapters::ResultDiagnostic{
                                                  .severity = adapters::DiagnosticSeverity::warning,
                                                  .code = import_diagnostic.code,
                                                  .message = import_diagnostic.message,
                                                  .molecule_index = molecule_index,
                                                  .line = import_diagnostic.line});
        }
        if (diagnostic.has_value()) {
            append_unique(record_diagnostics, *diagnostic);
        }
        if (result.status == calculation::ExecutionStatus::no_executable_plan) {
            for (const auto& rejected : result.rejections) {
                for (const auto& issue_value : rejected.issues) {
                    auto candidate = "method '" + rejected.method_id + "'";
                    if (rejected.parameter_set_id.has_value()) {
                        candidate += ", parameter set '" + *rejected.parameter_set_id + "'";
                    }
                    if (const auto* issue = std::get_if<methods::PrerequisiteIssue>(&issue_value)) {
                        if (issue->molecule_index.has_value() &&
                            *issue->molecule_index != molecule_index) {
                            continue;
                        }
                        append_unique(record_diagnostics,
                                      adapters::ResultDiagnostic{
                                          .severity = adapters::DiagnosticSeverity::error,
                                          .code = std::string{prerequisite_code(issue->kind)},
                                          .message = candidate + ": " + issue->message,
                                          .molecule_index = issue->molecule_index,
                                          .atom_index = issue->atom_index,
                                          .bond_index = issue->bond_index,
                                          .conformer_index = issue->conformer_index});
                    } else {
                        const auto& execution_issue =
                            std::get<methods::ExecutionIssue>(issue_value);
                        if (execution_issue.molecule_index.has_value() &&
                            *execution_issue.molecule_index != molecule_index) {
                            continue;
                        }
                        append_unique(
                            record_diagnostics,
                            adapters::ResultDiagnostic{
                                .severity = adapters::DiagnosticSeverity::error,
                                .code = std::string{methods::to_string(execution_issue.kind)},
                                .message = candidate + ": " + execution_issue.message,
                                .molecule_index = execution_issue.molecule_index});
                    }
                }
            }
        }
        document.records.push_back(
            adapters::ChargeResultRecord{.identity = record.identity,
                                         .charges = result.charges,
                                         .status = status,
                                         .diagnostics = std::move(record_diagnostics)});
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

auto make_requested_provenance(const ImportedExportContext& export_context,
                               const calculation::AssessmentRequest& request)
    -> adapters::RequestedCalculationProvenance {
    auto requested = adapters::RequestedCalculationProvenance{
        .method_id = request.method_id,
        .parameter_set_id = request.parameter_set_id,
        .permissive_types = request.classification_options.permissive_types,
        .cutoff_atom_threshold = request.resource_policy.cutoff_atom_threshold,
        .cover_atom_threshold = request.resource_policy.cover_atom_threshold,
        .max_threads = request.resource_policy.max_threads,
        .execution_kind = std::string{calculation::to_string(request.execution_selection.kind())},
        .execution_radius = request.execution_selection.radius(),
        .execution_charge_correction = request.execution_selection.charge_correction().transform(
            [](const calculation::ChargeCorrectionPolicy policy) {
                return std::string{calculation::to_string(policy)};
            }),
        .structural_input_policy =
            export_context.structural_input_policy.has_value()
                ? std::optional{adapters::StructuralInputPolicyProvenance{
                      .selection = export_context.structural_input_policy->selection,
                      .bonds = export_context.structural_input_policy->bonds}}
                : std::nullopt,
        .conformer_selection = export_context.conformer_selection};
    for (const auto& [method_id, options] : request.method_options) {
        requested.method_options.emplace(method_id, options);
    }
    return requested;
}

auto write_calculation_outputs(const std::string& output_directory, const std::string& input_path,
                               const ImportedExportContext& export_context,
                               const adapters::RequestedCalculationProvenance& requested,
                               const calculation::ExecutionResult& result, CalculationRun& run)
    -> int {
    run.metrics.peak_resident_memory_mb = peak_resident_memory_mb();
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
    if (!result.calculated()) {
        run.metrics.ended_at = utc_timestamp();
        run.metrics.runtime_seconds =
            std::chrono::duration<double>{std::chrono::steady_clock::now() - run.started}.count();
        const auto document = result_document(export_context, requested, result, run.metrics);
        write_json(prefix.string() + ".json", document);
        report_diagnostics(document);
        switch (result.status) {
        case calculation::ExecutionStatus::invalid_input_or_request:
            return 2;
        case calculation::ExecutionStatus::no_executable_plan:
            return 3;
        case calculation::ExecutionStatus::numerical_failure:
            return 4;
        case calculation::ExecutionStatus::cancelled:
            return 5;
        case calculation::ExecutionStatus::success:
            throw std::logic_error{"successful result has no charges"};
        }
    }

    const auto structural_output = export_context.format == ImportedExportContext::Format::pdb ||
                                   export_context.format == ImportedExportContext::Format::mmcif;
    if (!structural_output && export_context.format == ImportedExportContext::Format::json &&
        std::ranges::any_of(export_context.records,
                            [](const adapters::ImportedMoleculeRecord& record) {
                                const auto& molecule = record.molecule;
                                return molecule.conformer_count() > 1;
                            })) {
        throw std::runtime_error{
            "JSON input with multiple conformers cannot be written to SDF or MOL2"};
    }

    const auto* charges = result.charges ? std::addressof(*result.charges) : nullptr;
    if (charges == nullptr) {
        throw std::runtime_error{"calculation result is missing charges"};
    }
    const auto writing_started = std::chrono::steady_clock::now();
    write_mmcif(prefix.string() + ".cif", export_context, *charges);
    if (structural_output) {
        std::println("Wrote {} and {}", prefix.string() + ".json", prefix.string() + ".cif");
    } else {
        const auto assignments = assignments_by_molecule(*charges, export_context.records.size());
        write_sdf(prefix.string() + ".sdf", input_path, export_context, *charges);
        write_mol2(prefix.string() + ".mol2", input_path, export_context, assignments);
        std::println("Wrote {}, {}, {}, and {}", prefix.string() + ".json",
                     prefix.string() + ".sdf", prefix.string() + ".mol2", prefix.string() + ".cif");
    }
    run.metrics.writing_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - writing_started}.count();
    run.metrics.peak_resident_memory_mb = peak_resident_memory_mb();
    run.metrics.ended_at = utc_timestamp();
    run.metrics.runtime_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - run.started}.count();
    const auto document = result_document(export_context, requested, result, run.metrics);
    write_json(prefix.string() + ".json", document);
    report_diagnostics(document);
    return 0;
}

} // namespace chargefw::cli
