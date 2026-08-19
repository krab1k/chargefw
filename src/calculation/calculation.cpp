#include <chargefw/calculation/calculation.h>

#include "calculation/cutoff_execution.h"

#include <chargefw/methods/method.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_registry.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::calculation {
namespace {

[[nodiscard]] auto parameter_priority_of(const methods::ApplicableMethod& candidate) noexcept
    -> unsigned int {
    if (candidate.parameter_set == nullptr) {
        return 0;
    }

    return candidate.parameter_set->priority();
}

[[nodiscard]] auto parameter_id_of(const methods::ApplicableMethod& candidate) noexcept
    -> std::string_view {
    if (candidate.parameter_set == nullptr) {
        return {};
    }

    return candidate.parameter_set->id();
}

[[nodiscard]] auto ranks_before(const methods::ApplicableMethod& first,
                                const methods::ApplicableMethod& second) noexcept -> bool {
    const auto first_method_priority = first.method->metadata().priority;
    const auto second_method_priority = second.method->metadata().priority;

    if (first_method_priority != second_method_priority) {
        return first_method_priority > second_method_priority;
    }

    const auto first_parameter_priority = parameter_priority_of(first);
    const auto second_parameter_priority = parameter_priority_of(second);

    if (first_parameter_priority != second_parameter_priority) {
        return first_parameter_priority > second_parameter_priority;
    }

    if (first.method->id() != second.method->id()) {
        return first.method->id() < second.method->id();
    }

    return parameter_id_of(first) < parameter_id_of(second);
}

[[nodiscard]] auto assessment_for(const methods::ApplicableMethod& candidate,
                                  const ExecutionMode mode) -> const methods::ExecutionAssessment* {
    const auto found =
        std::ranges::find_if(candidate.execution_assessments,
                             [mode](const methods::ExecutionAssessment& assessment) -> bool {
                                 return assessment.mode == mode;
                             });
    return found == candidate.execution_assessments.end() ? nullptr : &*found;
}

[[nodiscard]] auto ranked_candidates(const methods::ApplicabilityResult& applicability)
    -> std::vector<const methods::ApplicableMethod*> {
    auto candidates = std::vector<const methods::ApplicableMethod*>{};
    candidates.reserve(applicability.applicable.size());

    for (const auto& candidate : applicability.applicable) {
        candidates.push_back(&candidate);
    }

    std::ranges::sort(candidates, [](const auto* first, const auto* second) -> bool {
        return ranks_before(*first, *second);
    });
    return candidates;
}

[[nodiscard]] auto plan_for(const methods::ApplicableMethod& candidate, const ExecutionMode mode,
                            const std::optional<double> radius,
                            const ChargeCorrectionPolicy charge_correction,
                            const methods::ExecutionAssessment& assessment) -> ExecutionPlan {
    return ExecutionPlan{.selected = &candidate,
                         .policy = ExecutionPolicy{mode, radius, charge_correction},
                         .issues = assessment.issues};
}

[[nodiscard]] auto application_methods(const ApplicationCalculationRequest& request)
    -> std::vector<const methods::Method*> {
    const auto& registry = methods::method_registry();

    if (request.method_id.has_value()) {
        const auto* method = registry.find(*request.method_id);

        if (method == nullptr) {
            throw std::invalid_argument{"method '" + *request.method_id + "' is not registered"};
        }

        return {method};
    }

    std::vector<const methods::Method*> result;
    result.reserve(registry.methods().size());

    for (const auto& method : registry.methods()) {
        result.push_back(method.get());
    }

    return result;
}

[[nodiscard]] auto application_parameter_sets(const ApplicationCalculationRequest& request)
    -> std::vector<parameters::ParameterSet> {
    if (!request.parameter_set_id.has_value()) {
        return request.parameter_sets;
    }

    const auto found = std::ranges::find_if(
        request.parameter_sets, [&request](const parameters::ParameterSet& parameter_set) -> bool {
            return parameter_set.id() == *request.parameter_set_id;
        });

    if (found == request.parameter_sets.end()) {
        throw std::invalid_argument{"parameter set '" + *request.parameter_set_id +
                                    "' was not provided"};
    }

    return {*found};
}

} // namespace

auto select_applicable_method(const methods::ApplicabilityResult& applicability)
    -> const methods::ApplicableMethod* {
    if (applicability.empty()) {
        return nullptr;
    }

    return &*std::ranges::min_element(applicability.applicable, ranks_before);
}

auto select_execution_plan(const methods::ApplicabilityResult& applicability,
                           const ExecutionSelection& selection) -> std::optional<ExecutionPlan> {
    const auto candidates = ranked_candidates(applicability);

    const auto select_mode =
        [&candidates, &selection](const ExecutionMode mode,
                                  const bool allow_warnings) -> std::optional<ExecutionPlan> {
        for (const auto* candidate : candidates) {
            const auto* assessment = assessment_for(*candidate, mode);
            if (assessment == nullptr ||
                assessment->availability == methods::ExecutionAvailability::unsupported ||
                (!allow_warnings && assessment->availability ==
                                        methods::ExecutionAvailability::available_with_warning)) {
                continue;
            }

            const auto radius =
                mode == ExecutionMode::full ? std::optional<double>{} : selection.radius();
            const auto charge_correction =
                mode == ExecutionMode::full
                    ? ChargeCorrectionPolicy::none
                    : selection.charge_correction().value_or(ChargeCorrectionPolicy::uniform);
            return plan_for(*candidate, mode, radius, charge_correction, *assessment);
        }

        return std::nullopt;
    };

    switch (selection.kind()) {
    case ExecutionSelectionKind::automatic:
        return select_mode(ExecutionMode::full, false);
    case ExecutionSelectionKind::full:
        return select_mode(ExecutionMode::full, true);
    case ExecutionSelectionKind::cutoff:
        return select_mode(ExecutionMode::cutoff, false);
    case ExecutionSelectionKind::cover:
        return select_mode(ExecutionMode::cover, false);
    }

    throw std::invalid_argument{"unknown execution selection"};
}

auto calculate(const CalculationRequest& request) -> CalculationResult {
    switch (request.execution_policy.mode()) {
    case ExecutionMode::full:
        return CalculationResult{
            .charges = methods::calculate_charges(request.selected, request.molecules)};
    case ExecutionMode::cutoff:
        return CalculationResult{.charges =
                                     calculate_cutoff_charges(request.selected, request.molecules,
                                                              request.execution_policy)};
    case ExecutionMode::cover:
        throw std::invalid_argument{"selected cover execution policy is not implemented"};
    }

    throw std::invalid_argument{"unknown execution policy"};
}

auto calculate(const ApplicationCalculationRequest& request) -> ApplicationCalculationResult {
    const auto candidate_methods = application_methods(request);
    const auto parameter_sets = application_parameter_sets(request);
    const features::PreparedMoleculeCollection prepared{request.molecules};

    auto applicability =
        methods::find_applicable_methods({.molecules = prepared,
                                          .methods = candidate_methods,
                                          .parameter_sets = parameter_sets,
                                          .classification_options = request.classification_options,
                                          .resource_policy = request.resource_policy});
    const auto plan = select_execution_plan(applicability, request.execution_selection);

    if (!plan.has_value() &&
        (request.method_id.has_value() || request.parameter_set_id.has_value() ||
         request.execution_selection.kind() != ExecutionSelectionKind::automatic)) {
        throw std::invalid_argument{"requested calculation selection has no executable plan"};
    }

    if (!plan.has_value()) {
        return ApplicationCalculationResult{.charges = std::nullopt,
                                            .applicability = std::move(applicability),
                                            .execution_policy = std::nullopt,
                                            .execution_issues = {}};
    }

    auto result = calculate(CalculationRequest{
        .molecules = prepared, .selected = *plan->selected, .execution_policy = plan->policy});
    return ApplicationCalculationResult{.charges = std::move(result.charges),
                                        .applicability = std::move(applicability),
                                        .execution_policy = plan->policy,
                                        .execution_issues = std::move(plan->issues)};
}

} // namespace chargefw::calculation
