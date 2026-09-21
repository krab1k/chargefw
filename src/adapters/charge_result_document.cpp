#include <chargefw/adapters/charge_result_document.h>

#include <chargefw/calculation/calculation.h>
#include <chargefw/methods/method_prerequisites.h>

#include <algorithm>
#include <functional>
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

auto validate_import_metadata(const std::span<const ImportedMoleculeRecord> records) -> void {
    for (const auto& record : records) {
        if (!record.import_metadata.has_value()) {
            continue;
        }
        const auto& metadata = *record.import_metadata;
        if (metadata.atoms.size() != record.molecule.atom_count()) {
            throw std::invalid_argument{
                "import atom mapping size does not match molecule atom count"};
        }
        if (metadata.conformers.size() != record.molecule.conformer_count()) {
            throw std::invalid_argument{
                "import conformer mapping size does not match molecule conformer count"};
        }
        for (const auto& conformer : metadata.conformers) {
            if (conformer.sites.size() != record.molecule.atom_count()) {
                throw std::invalid_argument{
                    "import conformer site mapping size does not match molecule atom count"};
            }
        }
    }
}

auto validate_assignments(const std::span<const ImportedMoleculeRecord> records,
                          const calculation::ExecutionResult& result) -> void {
    if (result.status != calculation::ExecutionStatus::success) {
        if (!result.charges.has_value()) {
            return;
        }
        throw std::invalid_argument{"unsuccessful calculation result contains charge assignments"};
    }
    if (!result.charges.has_value()) {
        throw std::invalid_argument{"successful calculation result has no charge assignments"};
    }
    if (!result.effective.has_value()) {
        throw std::invalid_argument{"successful calculation result has no effective provenance"};
    }
    if (result.charges->method_id() != result.effective->method_id) {
        throw std::invalid_argument{"charge assignment method does not match effective provenance"};
    }
    const auto charge_parameter_set = result.charges->parameter_set_id();
    if (charge_parameter_set.has_value() != result.effective->parameter_set_id.has_value() ||
        (charge_parameter_set.has_value() &&
         *charge_parameter_set != *result.effective->parameter_set_id)) {
        throw std::invalid_argument{
            "charge assignment parameter set does not match effective provenance"};
    }

    struct AssignmentCoverage {
        bool molecule = false;
        std::vector<bool> conformers;
    };

    auto coverage = std::vector<AssignmentCoverage>{};
    coverage.reserve(records.size());
    for (const auto& record : records) {
        coverage.push_back(
            AssignmentCoverage{.conformers = std::vector<bool>(record.molecule.conformer_count())});
    }

    auto conformer_scope = std::optional<bool>{};
    auto actual_order = std::vector<std::pair<std::size_t, std::optional<std::size_t>>>{};
    actual_order.reserve(result.charges->size());
    for (const auto& assignment : result.charges->assignments()) {
        const auto molecule_index = assignment.target.molecule_index;
        if (molecule_index >= records.size()) {
            throw std::invalid_argument{"charge assignment molecule index is outside the input"};
        }
        const auto& molecule = records[molecule_index].molecule;
        if (assignment.charges.size() != molecule.atom_count()) {
            throw std::invalid_argument{
                "charge assignment size does not match molecule atom count"};
        }

        const auto assignment_conformer_scope = assignment.target.conformer_index.has_value();
        if (conformer_scope.has_value() && *conformer_scope != assignment_conformer_scope) {
            throw std::invalid_argument{
                "calculation result mixes molecule and conformer assignment scopes"};
        }
        conformer_scope = assignment_conformer_scope;
        actual_order.emplace_back(molecule_index, assignment.target.conformer_index);

        auto& assigned = coverage[molecule_index];
        if (!assignment.target.conformer_index.has_value()) {
            if (assigned.molecule || std::ranges::any_of(assigned.conformers, std::identity{})) {
                throw std::invalid_argument{
                    "duplicate or mixed-scope charge assignments for one molecule"};
            }
            assigned.molecule = true;
            continue;
        }

        const auto conformer_index = *assignment.target.conformer_index;
        if (conformer_index >= molecule.conformer_count()) {
            throw std::invalid_argument{
                "charge assignment conformer index is outside the molecule"};
        }
        if (assigned.molecule || assigned.conformers[conformer_index]) {
            throw std::invalid_argument{
                "duplicate or mixed-scope charge assignments for one molecule"};
        }
        assigned.conformers[conformer_index] = true;
    }

    for (const auto& assigned : coverage) {
        if (!assigned.molecule && (assigned.conformers.empty() ||
                                   !std::ranges::all_of(assigned.conformers, std::identity{}))) {
            throw std::invalid_argument{"calculation result does not cover every input molecule"};
        }
    }

    auto expected_order = std::vector<std::pair<std::size_t, std::optional<std::size_t>>>{};
    if (conformer_scope.value_or(false)) {
        for (std::size_t molecule_index = 0; molecule_index < records.size(); ++molecule_index) {
            for (std::size_t conformer_index = 0;
                 conformer_index < records[molecule_index].molecule.conformer_count();
                 ++conformer_index) {
                expected_order.emplace_back(molecule_index, conformer_index);
            }
        }
    } else {
        for (std::size_t molecule_index = 0; molecule_index < records.size(); ++molecule_index) {
            expected_order.emplace_back(molecule_index, std::nullopt);
        }
    }
    if (actual_order != expected_order) {
        throw std::invalid_argument{"charge assignments are not in canonical input order"};
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

auto make_charge_calculation_result(std::vector<ImportedMoleculeRecord> inputs,
                                    RequestedCalculationProvenance requested,
                                    calculation::ExecutionResult execution)
    -> ChargeCalculationResult {
    validate_import_metadata(inputs);
    validate_assignments(inputs, execution);
    auto result = ChargeCalculationResult{};
    result.inputs_ = std::move(inputs);
    result.requested_ = std::move(requested);
    result.execution_ = std::move(execution);
    return result;
}

auto ChargeCalculationResult::inputs() const noexcept -> std::span<const ImportedMoleculeRecord> {
    return inputs_;
}

auto ChargeCalculationResult::requested() const noexcept -> const RequestedCalculationProvenance& {
    return requested_;
}

auto ChargeCalculationResult::execution() const noexcept -> const calculation::ExecutionResult& {
    return execution_;
}

auto make_charge_result_document(const ChargeCalculationResult& result,
                                 const std::string_view generator_name,
                                 const std::string_view generator_version,
                                 std::optional<ExecutionMetrics> execution_metrics)
    -> ChargeResultDocument {
    const auto& execution = result.execution();
    const auto records = result.inputs();
    const auto diagnostic = calculation_diagnostic(execution);
    auto effective = EffectiveCalculationProvenance{};
    if (execution.effective.has_value()) {
        const auto& calculation = *execution.effective;
        effective.method_id = calculation.method_id;
        effective.parameter_set_id = calculation.parameter_set_id;
        effective.execution_mode =
            std::string{calculation::to_string(calculation.execution_policy.mode())};
        effective.execution_radius = calculation.execution_policy.radius();
        effective.warnings.reserve(calculation.execution_issues.size());
        for (const auto& issue : calculation.execution_issues) {
            effective.warnings.push_back(issue.message);
        }
        effective.method_options.emplace(calculation.method_id, calculation.method_options);
    }
    auto document = ChargeResultDocument{
        .generator_name = std::string{generator_name},
        .generator_version = std::string{generator_version},
        .status = execution.status,
        .diagnostics =
            diagnostic.has_value() ? std::vector{*diagnostic} : std::vector<ResultDiagnostic>{},
        .records = {},
        .calculation_provenance =
            CalculationProvenance{.requested = result.requested(),
                                  .effective = std::move(effective),
                                  .execution_metrics = std::move(execution_metrics)}};
    document.records.reserve(records.size());
    for (std::size_t molecule_index = 0; molecule_index < records.size(); ++molecule_index) {
        const auto& record = records[molecule_index];
        auto assignments = std::vector<charges::ChargeAssignment>{};
        if (execution.charges.has_value()) {
            for (const auto& assignment : execution.charges->assignments()) {
                if (assignment.target.molecule_index == molecule_index) {
                    assignments.push_back(assignment);
                }
            }
        }
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
        if (execution.status == calculation::ExecutionStatus::no_executable_plan) {
            for (const auto& rejected : execution.rejections) {
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
            ChargeResultRecord{.input = record,
                               .assignments = std::move(assignments),
                               .status = execution.status,
                               .diagnostics = std::move(record_diagnostics)});
    }
    return document;
}

} // namespace chargefw::adapters
