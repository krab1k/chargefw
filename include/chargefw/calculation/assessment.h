#pragma once

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/parameters/classification/classification_result.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <memory>
#include <optional>
#include <vector>

namespace chargefw::calculation {

struct ApplicationCalculationRequest;

// Non-owning concrete execution choice produced from an applicability result. The selected
// candidate and its classifications remain owned by the ApplicabilityResult.
struct ExecutionPlan {
    const methods::ApplicableMethod* selected = nullptr;
    ExecutionPolicy policy{};
    std::vector<methods::ExecutionIssue> issues;
};

// Owns prepared application inputs and their applicability/execution-plan assessment. The contained
// candidates reference methods from the registry and the owned parameter sets below. Move this
// result into calculate() to execute without repeating preparation or classification.
class ApplicationAssessmentResult {
  public:
    ApplicationAssessmentResult(core::MoleculeCollection molecules,
                                std::vector<parameters::ParameterSet> parameter_sets);
    ApplicationAssessmentResult(const ApplicationAssessmentResult&) = delete;
    auto operator=(const ApplicationAssessmentResult&) -> ApplicationAssessmentResult& = delete;
    ApplicationAssessmentResult(ApplicationAssessmentResult&&) noexcept = default;
    auto operator=(ApplicationAssessmentResult&&) noexcept
        -> ApplicationAssessmentResult& = default;
    ~ApplicationAssessmentResult();

    [[nodiscard]] auto prepared_molecules() const noexcept
        -> const features::PreparedMoleculeCollection&;

    std::vector<parameters::ParameterSet> parameter_sets;
    methods::ApplicabilityResult applicability;
    const methods::ApplicableMethod* selected = nullptr;
    std::optional<ExecutionPolicy> execution_policy;
    std::vector<methods::ExecutionIssue> execution_issues;
    double applicability_seconds = 0.0;

    [[nodiscard]] auto executable() const noexcept -> bool {
        return execution_policy.has_value();
    }

  private:
    std::unique_ptr<core::MoleculeCollection> molecules_;
    std::unique_ptr<features::PreparedMoleculeCollection> prepared_molecules_;
};

// Selects the applicable candidate with the highest method priority, then the highest parameter-set
// priority. Equal priorities are resolved by method ID, then parameter-set ID, in lexicographic
// order. Returns nullptr if no candidate is applicable.
[[nodiscard]] auto select_applicable_method(const methods::ApplicabilityResult& applicability)
    -> const methods::ApplicableMethod*;

// Selects a concrete execution plan from scientifically applicable candidates. Automatic selection
// considers only non-discouraged full execution until a reduced executor is implemented.
[[nodiscard]] auto select_execution_plan(const methods::ApplicabilityResult& applicability,
                                         const ExecutionSelection& selection)
    -> std::optional<ExecutionPlan>;

// Resolves registered methods, applies classification policy, and selects a concrete execution plan
// without running a charge calculation. Explicit unavailable method or parameter-set IDs are
// errors; an unavailable explicit execution plan is reported as no plan.
[[nodiscard]] auto assess(const ApplicationCalculationRequest& request)
    -> ApplicationAssessmentResult;

} // namespace chargefw::calculation
