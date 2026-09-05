#include <chargefw/adapters/charge_result_document.h>

#include <chargefw/calculation/calculation.h>
#include <chargefw/methods/method_prerequisites.h>

#include <algorithm>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace chargefw::adapters {
namespace {

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

auto append_unique(std::vector<ResultDiagnostic>& diagnostics, ResultDiagnostic diagnostic)
    -> void {
    const auto duplicate =
        std::ranges::any_of(diagnostics, [&diagnostic](const ResultDiagnostic& existing) {
            return existing.code == diagnostic.code && existing.message == diagnostic.message;
        });
    if (!duplicate) {
        diagnostics.push_back(std::move(diagnostic));
    }
}

[[nodiscard]] auto calculation_diagnostic(const calculation::ExecutionResult& result)
    -> std::optional<ResultDiagnostic> {
    const auto make = [](const DiagnosticSeverity severity, std::string code, std::string message) {
        auto diagnostic = ResultDiagnostic{};
        diagnostic.severity = severity;
        diagnostic.code = std::move(code);
        diagnostic.message = std::move(message);
        return diagnostic;
    };
    switch (result.status) {
    case calculation::ExecutionStatus::success:
        return std::nullopt;
    case calculation::ExecutionStatus::invalid_input_or_request:
        return make(DiagnosticSeverity::error, "invalid_input_or_request",
                    result.failure_message.value_or("Invalid input or calculation request."));
    case calculation::ExecutionStatus::no_executable_plan:
        return make(DiagnosticSeverity::error, "no_executable_plan",
                    "No executable method and parameter-set plan was found.");
    case calculation::ExecutionStatus::numerical_failure:
        return make(DiagnosticSeverity::error, "numerical_failure",
                    result.failure_message.value_or("Calculation failed numerically."));
    case calculation::ExecutionStatus::cancelled:
        return make(DiagnosticSeverity::info, "calculation_cancelled",
                    "Calculation was cancelled before completion.");
    }
    throw std::logic_error{"unknown execution result status"};
}

} // namespace

auto make_charge_result_document(const std::span<const ImportedMoleculeRecord> records,
                                 const RequestedCalculationProvenance& requested,
                                 const calculation::ExecutionResult& result,
                                 const std::string_view generator_name,
                                 const std::string_view generator_version,
                                 std::optional<ExecutionMetrics> execution_metrics)
    -> ChargeResultDocument {
    const auto diagnostic = calculation_diagnostic(result);
    auto effective = EffectiveCalculationProvenance{};
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
    auto document = ChargeResultDocument{
        .generator_name = std::string{generator_name},
        .generator_version = std::string{generator_version},
        .status = result.status,
        .diagnostics =
            diagnostic.has_value() ? std::vector{*diagnostic} : std::vector<ResultDiagnostic>{},
        .records = {},
        .calculation_provenance =
            CalculationProvenance{.requested = requested,
                                  .effective = std::move(effective),
                                  .execution_metrics = std::move(execution_metrics)}};
    document.records.reserve(records.size());
    for (std::size_t molecule_index = 0; molecule_index < records.size(); ++molecule_index) {
        const auto& record = records[molecule_index];
        auto record_diagnostics = std::vector<ResultDiagnostic>{};
        for (const auto& import_diagnostic : record.diagnostics) {
            auto imported = ResultDiagnostic{};
            imported.severity = DiagnosticSeverity::warning;
            imported.code = import_diagnostic.code;
            imported.message = import_diagnostic.message;
            imported.molecule_index = molecule_index;
            imported.line = import_diagnostic.line;
            append_unique(record_diagnostics, std::move(imported));
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
                        auto rejected_issue = ResultDiagnostic{};
                        rejected_issue.severity = DiagnosticSeverity::error;
                        rejected_issue.code = prerequisite_code(issue->kind);
                        rejected_issue.message = candidate + ": " + issue->message;
                        rejected_issue.molecule_index = issue->molecule_index;
                        rejected_issue.atom_index = issue->atom_index;
                        rejected_issue.bond_index = issue->bond_index;
                        rejected_issue.conformer_index = issue->conformer_index;
                        append_unique(record_diagnostics, std::move(rejected_issue));
                    } else {
                        const auto& execution_issue =
                            std::get<methods::ExecutionIssue>(issue_value);
                        if (execution_issue.molecule_index.has_value() &&
                            *execution_issue.molecule_index != molecule_index) {
                            continue;
                        }
                        auto rejected_issue = ResultDiagnostic{};
                        rejected_issue.severity = DiagnosticSeverity::error;
                        rejected_issue.code = methods::to_string(execution_issue.kind);
                        rejected_issue.message = candidate + ": " + execution_issue.message;
                        rejected_issue.molecule_index = execution_issue.molecule_index;
                        append_unique(record_diagnostics, std::move(rejected_issue));
                    }
                }
            }
        }
        document.records.push_back(
            ChargeResultRecord{.identity = record.identity,
                               .charges = result.charges,
                               .status = result.status,
                               .diagnostics = std::move(record_diagnostics)});
    }
    return document;
}

} // namespace chargefw::adapters
