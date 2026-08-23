#pragma once

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/parameters/classification/classification_result.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

namespace chargefw::calculation {

class AssessmentResult;
class CalculationObserver;
struct AssessmentRequest;
struct ExecutionResult;
// Copies the owned molecule and parameter inputs into the assessment result.
[[nodiscard]] auto assess(const AssessmentRequest& request) -> AssessmentResult;
// Snapshots selection configuration and transfers execution inputs into the assessment result.
// Callers must not inspect the request after passing it to this consuming overload.
[[nodiscard]] auto assess(AssessmentRequest&& request) -> AssessmentResult;

// Non-owning concrete execution choice produced from an applicability result. The selected
// candidate and its classifications remain owned by the ApplicabilityResult.
struct ExecutionPlan {
    const methods::ApplicableMethod* selected = nullptr;
    ExecutionPolicy policy{};
    std::vector<methods::ExecutionIssue> issues;
};

// Owns application assessment inputs so adapters and bindings do not need to manage native method
// pointers, parameter spans, or prepared-feature lifetimes.
struct AssessmentRequest {
    core::MoleculeCollection molecules;
    std::vector<parameters::ParameterSet> parameter_sets;
    std::optional<std::string> method_id;
    std::optional<std::string> parameter_set_id;
    std::unordered_map<std::string, methods::MethodOptions> method_options;
    parameters::ClassificationOptions classification_options{};
    ExecutionSelection execution_selection{};
    ResourcePolicy resource_policy{};
};

// Value-only applicability data for application callers. It deliberately contains no registry or
// parameter-set pointers, so it remains valid after an assessment is consumed by calculate().
struct ApplicableCandidateReport {
    std::string method_id;
    std::optional<std::string> parameter_set_id;
    std::vector<methods::ExecutionAssessment> execution_assessments;
};

struct RejectedCandidateReport {
    std::string method_id;
    std::optional<std::string> parameter_set_id;
    std::vector<methods::PrerequisiteIssue> issues;
};

struct ApplicabilityReport {
    std::vector<ApplicableCandidateReport> applicable;
    std::vector<RejectedCandidateReport> rejected;
    std::optional<std::size_t> selected_candidate_index;

    [[nodiscard]] auto empty() const noexcept -> bool {
        return applicable.empty();
    }
};

// Owns prepared application inputs and their applicability/execution-plan assessment. Move this
// result into calculate() to execute without repeating preparation or classification.
class AssessmentResult {
  public:
    AssessmentResult(core::MoleculeCollection molecules,
                     std::vector<parameters::ParameterSet> parameter_sets,
                     bool requires_executable_plan);
    AssessmentResult(const AssessmentResult&) = delete;
    auto operator=(const AssessmentResult&) -> AssessmentResult& = delete;
    AssessmentResult(AssessmentResult&&) noexcept = default;
    // Prepared features non-owningly refer to molecules_, so replacing the two independent owners
    // cannot preserve their lifetime invariant.
    auto operator=(AssessmentResult&&) noexcept -> AssessmentResult& = delete;
    ~AssessmentResult();

    [[nodiscard]] auto prepared_molecules() const noexcept
        -> const features::PreparedMoleculeCollection&;

    [[nodiscard]] auto applicability() const noexcept -> const ApplicabilityReport&;
    [[nodiscard]] auto execution_policy() const noexcept -> const std::optional<ExecutionPolicy>&;
    [[nodiscard]] auto execution_issues() const noexcept
        -> const std::vector<methods::ExecutionIssue>&;
    [[nodiscard]] auto applicability_seconds() const noexcept -> double;

    [[nodiscard]] auto executable() const noexcept -> bool {
        return execution_policy_.has_value();
    }

    [[nodiscard]] auto requires_executable_plan() const noexcept -> bool {
        return requires_executable_plan_;
    }

  private:
    auto
    assess_prepared(std::span<const methods::Method* const> selected_methods,
                    const parameters::ClassificationOptions& classification_options,
                    const ResourcePolicy& resource_policy,
                    const std::unordered_map<std::string, methods::MethodOptions>& method_options,
                    const ExecutionSelection& execution_selection) -> void;

    friend auto assess(const AssessmentRequest& request) -> AssessmentResult;
    friend auto assess(AssessmentRequest&& request) -> AssessmentResult;
    friend auto calculate(AssessmentResult assessment, std::size_t max_threads,
                          const CalculationObserver& observer) -> ExecutionResult;

    std::vector<parameters::ParameterSet> parameter_sets_;
    methods::ApplicabilityResult applicability_;
    ApplicabilityReport applicability_report_;
    std::optional<std::size_t> selected_candidate_index_;
    std::optional<ExecutionPolicy> execution_policy_;
    std::vector<methods::ExecutionIssue> execution_issues_;
    double applicability_seconds_ = 0.0;
    bool requires_executable_plan_ = false;
    std::unique_ptr<core::MoleculeCollection> molecules_;
    std::unique_ptr<features::PreparedMoleculeCollection> prepared_molecules_;
};

// Selects the applicable candidate with the highest method priority, then the highest parameter-set
// priority. Equal priorities are resolved by method ID, then parameter-set ID, in lexicographic
// order. Returns nullptr if no candidate is applicable.
[[nodiscard]] auto select_applicable_method(const methods::ApplicabilityResult& applicability)
    -> const methods::ApplicableMethod*;

// Selects a concrete execution plan from scientifically applicable candidates.
[[nodiscard]] auto select_execution_plan(const methods::ApplicabilityResult& applicability,
                                         const ExecutionSelection& selection)
    -> std::optional<ExecutionPlan>;

} // namespace chargefw::calculation
