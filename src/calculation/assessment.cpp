#include <chargefw/calculation/calculation.h>

#include <chargefw/methods/method.h>
#include <chargefw/methods/method_registry.h>

#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::calculation {
namespace {

[[nodiscard]] auto parameter_priority_of(const methods::ApplicableMethod& candidate) noexcept
    -> unsigned int {
    return candidate.parameter_set == nullptr ? 0U : candidate.parameter_set->priority();
}

[[nodiscard]] auto parameter_id_of(const methods::ApplicableMethod& candidate) noexcept
    -> std::string_view {
    return candidate.parameter_set == nullptr ? std::string_view{} : candidate.parameter_set->id();
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
    const auto found = std::ranges::find_if(
        candidate.execution_assessments,
        [mode](const methods::ExecutionAssessment& assessment) { return assessment.mode == mode; });
    return found == candidate.execution_assessments.end() ? nullptr : &*found;
}

[[nodiscard]] auto ranked_candidates(const methods::ApplicabilityResult& applicability)
    -> std::vector<const methods::ApplicableMethod*> {
    auto candidates = std::vector<const methods::ApplicableMethod*>{};
    candidates.reserve(applicability.applicable.size());
    for (const auto& candidate : applicability.applicable) {
        candidates.push_back(&candidate);
    }
    std::ranges::sort(candidates, [](const auto* first, const auto* second) {
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

[[nodiscard]] auto application_methods(const AssessmentRequest& request)
    -> std::vector<const methods::Method*> {
    const auto& registry = methods::method_registry();
    if (request.method_id.has_value()) {
        const auto* method = registry.find(*request.method_id);
        if (method == nullptr) {
            throw std::invalid_argument{"method '" + *request.method_id + "' is not registered"};
        }
        return {method};
    }

    auto result = std::vector<const methods::Method*>{};
    result.reserve(registry.methods().size());
    for (const auto& method : registry.methods()) {
        result.push_back(method.get());
    }
    return result;
}

auto validate_application_method_options(const AssessmentRequest& request) -> void {
    const auto& registry = methods::method_registry();
    for (const auto& [method_id, overrides] : request.method_options) {
        const auto* method = registry.find(method_id);
        if (method == nullptr) {
            throw std::invalid_argument{"method '" + method_id + "' is not registered"};
        }
        if (request.method_id.has_value() && method_id != *request.method_id) {
            throw std::invalid_argument{"method options supplied for '" + method_id +
                                        "' but method '" + *request.method_id + "' was requested"};
        }

        auto options = methods::make_default_options(method->option_schema());
        for (const auto& [id, value] : overrides.values()) {
            options.set(id, value);
        }
        methods::validate_method_options(method->option_schema(), options);
    }
}

[[nodiscard]] auto application_parameter_sets(const AssessmentRequest& request)
    -> std::vector<parameters::ParameterSet> {
    if (!request.parameter_set_id.has_value()) {
        return request.parameter_sets;
    }

    const auto found = std::ranges::find_if(
        request.parameter_sets, [&request](const parameters::ParameterSet& parameter_set) {
            return parameter_set.id() == *request.parameter_set_id;
        });
    if (found == request.parameter_sets.end()) {
        throw std::invalid_argument{"parameter set '" + *request.parameter_set_id +
                                    "' was not provided"};
    }
    return {*found};
}

auto assess_prepared(const AssessmentRequest& request, AssessmentResult& result) -> void {
    result.applicability =
        methods::find_applicable_methods({.molecules = result.prepared_molecules(),
                                          .methods = application_methods(request),
                                          .parameter_sets = result.parameter_sets,
                                          .classification_options = request.classification_options,
                                          .resource_policy = request.resource_policy,
                                          .method_options = request.method_options});
    const auto plan = select_execution_plan(result.applicability, request.execution_selection);
    if (plan.has_value()) {
        result.execution_policy = plan->policy;
        result.selected = plan->selected;
        result.execution_issues = plan->issues;
    }
}

} // namespace

AssessmentResult::AssessmentResult(core::MoleculeCollection molecules,
                                   std::vector<parameters::ParameterSet> supplied_parameter_sets,
                                   const bool requires_executable_plan)
    : parameter_sets{std::move(supplied_parameter_sets)},
      requires_executable_plan_{requires_executable_plan},
      molecules_{std::make_unique<core::MoleculeCollection>(std::move(molecules))},
      prepared_molecules_{std::make_unique<features::PreparedMoleculeCollection>(*molecules_)} {}

AssessmentResult::~AssessmentResult() = default;

auto AssessmentResult::prepared_molecules() const noexcept
    -> const features::PreparedMoleculeCollection& {
    return *prepared_molecules_;
}

auto select_applicable_method(const methods::ApplicabilityResult& applicability)
    -> const methods::ApplicableMethod* {
    return applicability.empty()
               ? nullptr
               : &*std::ranges::min_element(applicability.applicable, ranks_before);
}

auto select_execution_plan(const methods::ApplicabilityResult& applicability,
                           const ExecutionSelection& selection) -> std::optional<ExecutionPlan> {
    const auto candidates = ranked_candidates(applicability);
    const auto select_mode = [&candidates, &selection](const ExecutionMode mode,
                                                       const bool allow_warnings) {
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
            const auto correction =
                mode == ExecutionMode::full
                    ? ChargeCorrectionPolicy::none
                    : selection.charge_correction().value_or(ChargeCorrectionPolicy::uniform);
            return std::optional{plan_for(*candidate, mode, radius, correction, *assessment)};
        }
        return std::optional<ExecutionPlan>{};
    };

    switch (selection.kind()) {
    case ExecutionSelectionKind::automatic: {
        const auto radius = selection.radius().value_or(default_automatic_reduced_radius);
        for (const auto* candidate : candidates) {
            const auto* full = assessment_for(*candidate, ExecutionMode::full);
            if (full != nullptr &&
                full->availability == methods::ExecutionAvailability::available) {
                return plan_for(*candidate, ExecutionMode::full, std::nullopt,
                                ChargeCorrectionPolicy::none, *full);
            }
            const auto* cutoff = assessment_for(*candidate, ExecutionMode::cutoff);
            if (cutoff != nullptr &&
                cutoff->availability != methods::ExecutionAvailability::unsupported) {
                return plan_for(*candidate, ExecutionMode::cutoff, radius,
                                ChargeCorrectionPolicy::uniform, *cutoff);
            }
            const auto* cover = assessment_for(*candidate, ExecutionMode::cover);
            if (cover != nullptr &&
                cover->availability != methods::ExecutionAvailability::unsupported) {
                return plan_for(*candidate, ExecutionMode::cover, radius,
                                ChargeCorrectionPolicy::uniform, *cover);
            }
        }
        return std::nullopt;
    }
    case ExecutionSelectionKind::full:
        return select_mode(ExecutionMode::full, true);
    case ExecutionSelectionKind::cutoff:
        return select_mode(ExecutionMode::cutoff, false);
    case ExecutionSelectionKind::cover:
        return select_mode(ExecutionMode::cover, false);
    }
    throw std::invalid_argument{"unknown execution selection"};
}

auto assess(const AssessmentRequest& request) -> AssessmentResult {
    const auto started = std::chrono::steady_clock::now();
    validate_application_method_options(request);
    auto result = AssessmentResult{
        request.molecules, application_parameter_sets(request),
        request.method_id.has_value() || request.parameter_set_id.has_value() ||
            request.execution_selection.kind() != ExecutionSelectionKind::automatic};
    assess_prepared(request, result);
    result.applicability_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - started}.count();
    return result;
}

} // namespace chargefw::calculation
