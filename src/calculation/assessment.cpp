#include <chargefw/calculation/assessment.h>

#include <chargefw/methods/method.h>
#include <chargefw/methods/method_registry.h>

#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace chargefw::calculation {

class PlanIdentity {};

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

[[nodiscard]] auto assessment_methods(const AssessmentRequest& request)
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

auto validate_assessment_method_options(const AssessmentRequest& request) -> void {
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

auto validate_unique_parameter_set_ids(const AssessmentRequest& request) -> void {
    auto ids = std::unordered_set<std::string_view>{};
    ids.reserve(request.parameter_sets.size());
    for (const auto& parameter_set : request.parameter_sets) {
        const auto id = parameter_set.id();
        if (!ids.insert(id).second) {
            throw std::invalid_argument{"duplicate parameter set id '" + std::string{id} +
                                        "' in assessment request"};
        }
    }
}

auto validate_resource_policy(const ResourcePolicy& resource_policy) -> void {
    if (resource_policy.cover_atom_threshold.has_value() &&
        !resource_policy.cutoff_atom_threshold.has_value()) {
        throw std::invalid_argument{"cover atom threshold requires a finite cutoff atom threshold"};
    }
    if (resource_policy.cover_atom_threshold.has_value() &&
        *resource_policy.cover_atom_threshold < *resource_policy.cutoff_atom_threshold) {
        throw std::invalid_argument{
            "cover atom threshold must not be smaller than cutoff atom threshold"};
    }
}

auto retain_requested_parameter_set(AssessmentRequest& request) -> void {
    if (!request.parameter_set_id.has_value()) {
        return;
    }

    const auto found = std::ranges::find_if(
        request.parameter_sets, [&request](const parameters::ParameterSet& parameter_set) {
            return parameter_set.id() == *request.parameter_set_id;
        });
    if (found == request.parameter_sets.end()) {
        throw std::invalid_argument{"parameter set '" + *request.parameter_set_id +
                                    "' was not provided"};
    }

    std::erase_if(request.parameter_sets,
                  [&request](const parameters::ParameterSet& parameter_set) {
                      return parameter_set.id() != *request.parameter_set_id;
                  });
}

} // namespace

AssessmentResult::AssessmentResult(core::MoleculeCollection molecules,
                                   std::vector<parameters::ParameterSet> supplied_parameter_sets)
    : parameter_sets_{std::move(supplied_parameter_sets)},
      plan_identity_{std::make_shared<PlanIdentity>()},
      molecules_{std::make_unique<core::MoleculeCollection>(std::move(molecules))},
      prepared_molecules_{std::make_unique<features::PreparedMoleculeCollection>(*molecules_)} {}

AssessmentResult::~AssessmentResult() = default;

ExecutionPlan::ExecutionPlan(std::shared_ptr<const PlanIdentity> identity,
                             const methods::ApplicableMethod& candidate, ExecutionPolicy policy,
                             std::vector<methods::ExecutionIssue> warnings)
    : identity_{std::move(identity)}, candidate_{&candidate}, policy_{policy},
      warnings_{std::move(warnings)} {}

auto ExecutionPlan::candidate() const noexcept -> const methods::ApplicableMethod& {
    return *candidate_;
}

auto ExecutionPlan::policy() const noexcept -> const ExecutionPolicy& {
    return policy_;
}

auto ExecutionPlan::warnings() const noexcept -> std::span<const methods::ExecutionIssue> {
    return warnings_;
}

auto AssessmentResult::assess_prepared(
    const std::span<const methods::Method* const> selected_methods,
    const parameters::ClassificationOptions& classification_options,
    const ResourcePolicy& resource_policy,
    const std::unordered_map<std::string, methods::MethodOptions>& method_options,
    const ExecutionSelection& execution_selection) -> void {
    applicability_ =
        methods::find_applicable_methods({.molecules = prepared_molecules(),
                                          .methods = selected_methods,
                                          .parameter_sets = parameter_sets_,
                                          .classification_options = classification_options,
                                          .resource_policy = resource_policy,
                                          .method_options = method_options});
    rejections_.reserve(applicability_.rejected.size());
    for (const auto& candidate : applicability_.rejected) {
        auto issues = std::vector<RejectionIssue>{};
        issues.reserve(candidate.issues.size());
        for (const auto& issue : candidate.issues) {
            issues.emplace_back(issue);
        }
        rejections_.push_back(Rejection{
            .method_id = std::string{selected_methods[candidate.method_index]->id()},
            .parameter_set_id = candidate.parameter_set_index.has_value()
                                    ? std::optional{std::string{
                                          parameter_sets_[*candidate.parameter_set_index].id()}}
                                    : std::nullopt,
            .policy = std::nullopt,
            .issues = std::move(issues),
        });
    }
    const auto candidates = ranked_candidates(applicability_);
    const auto policy_for = [&execution_selection](const ExecutionMode mode) {
        const auto radius = mode == ExecutionMode::full
                                ? std::optional<double>{}
                                : std::optional{execution_selection.radius().value_or(
                                      default_automatic_reduced_radius)};
        const auto correction =
            mode == ExecutionMode::full
                ? ChargeCorrectionPolicy::none
                : execution_selection.charge_correction().value_or(ChargeCorrectionPolicy::uniform);
        return ExecutionPolicy{mode, radius, correction};
    };
    const auto append_rejection =
        [this, &policy_for](const methods::ApplicableMethod& candidate, const ExecutionMode mode,
                            const methods::ExecutionAssessment& assessment) {
            auto issues = std::vector<RejectionIssue>{};
            issues.reserve(assessment.issues.size());
            for (const auto& issue : assessment.issues) {
                issues.emplace_back(issue);
            }
            rejections_.push_back(Rejection{
                .method_id = std::string{candidate.method->id()},
                .parameter_set_id = candidate.parameter_set == nullptr
                                        ? std::nullopt
                                        : std::optional{std::string{candidate.parameter_set->id()}},
                .policy = policy_for(mode),
                .issues = std::move(issues),
            });
        };
    const auto consider_mode = [this, &policy_for, &append_rejection](
                                   const methods::ApplicableMethod& candidate,
                                   const ExecutionMode mode, const bool permit_warnings) {
        const auto* assessment = assessment_for(candidate, mode);
        if (assessment == nullptr) {
            return;
        }
        if (assessment->availability == methods::ExecutionAvailability::unsupported ||
            (!permit_warnings &&
             assessment->availability == methods::ExecutionAvailability::available_with_warning)) {
            append_rejection(candidate, mode, *assessment);
            return;
        }
        plans_.push_back(
            ExecutionPlan{plan_identity_, candidate, policy_for(mode), assessment->issues});
    };

    for (const auto* candidate : candidates) {
        switch (execution_selection.kind()) {
        case ExecutionSelectionKind::automatic:
            consider_mode(*candidate, ExecutionMode::full, false);
            consider_mode(*candidate, ExecutionMode::cutoff, false);
            consider_mode(*candidate, ExecutionMode::cover, false);
            break;
        case ExecutionSelectionKind::full:
            consider_mode(*candidate, ExecutionMode::full, true);
            break;
        case ExecutionSelectionKind::cutoff:
            consider_mode(*candidate, ExecutionMode::cutoff, true);
            break;
        case ExecutionSelectionKind::cover:
            consider_mode(*candidate, ExecutionMode::cover, true);
            break;
        }
    }
}

auto AssessmentResult::prepared_molecules() const noexcept
    -> const features::PreparedMoleculeCollection& {
    return *prepared_molecules_;
}

auto AssessmentResult::applicability_seconds() const noexcept -> double {
    return applicability_seconds_;
}

auto AssessmentResult::plans() const noexcept -> std::span<const ExecutionPlan> {
    return plans_;
}

auto AssessmentResult::rejections() const noexcept -> std::span<const Rejection> {
    return rejections_;
}

auto AssessmentResult::default_plan() const noexcept -> const ExecutionPlan* {
    return plans_.empty() ? nullptr : &plans_.front();
}

auto select_applicable_method(const methods::ApplicabilityResult& applicability)
    -> const methods::ApplicableMethod* {
    return applicability.empty()
               ? nullptr
               : &*std::ranges::min_element(applicability.applicable, ranks_before);
}

auto AssessmentResult::assess_owned(AssessmentRequest request) -> AssessmentResult {
    const auto started = std::chrono::steady_clock::now();
    validate_assessment_method_options(request);
    validate_unique_parameter_set_ids(request);
    validate_resource_policy(request.resource_policy);
    const auto selected_methods = assessment_methods(request);
    const auto classification_options = request.classification_options;
    const auto execution_selection = request.execution_selection;
    const auto resource_policy = request.resource_policy;
    retain_requested_parameter_set(request);
    auto method_options = std::move(request.method_options);
    auto result = AssessmentResult{std::move(request.molecules), std::move(request.parameter_sets)};
    result.assess_prepared(selected_methods, classification_options, resource_policy,
                           method_options, execution_selection);
    result.applicability_seconds_ =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - started}.count();
    return result;
}

auto assess(const AssessmentRequest& request) -> AssessmentResult {
    return AssessmentResult::assess_owned(request);
}

auto assess(AssessmentRequest&& request) -> AssessmentResult {
    return AssessmentResult::assess_owned(std::move(request));
}

} // namespace chargefw::calculation
