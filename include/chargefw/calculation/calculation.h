#pragma once

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/parameters/classification/classification_result.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace chargefw::calculation {

struct CalculationRequest {
    // Non-owning input; the referenced collection and selected candidate must outlive calculate().
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
    const features::PreparedMoleculeCollection& molecules;
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
    const methods::ApplicableMethod& selected;
    ExecutionPolicy execution_policy{};
    // Zero delegates the worker count to the oneTBB runtime.
    std::size_t max_threads = 0;
    // Optional, non-owning observer for progress and cancellation. Nullptr means no observation.
    const CalculationObserver* observer = nullptr;
};

struct CalculationResult {
    charges::ChargeSet charges;
};

struct CalculationMetrics {
    double applicability_seconds = 0.0;
    double computation_seconds = 0.0;
};

// Non-owning concrete execution choice produced from an applicability result. The selected
// candidate and its classifications remain owned by the ApplicabilityResult.
struct ExecutionPlan {
    const methods::ApplicableMethod* selected = nullptr;
    ExecutionPolicy policy{};
    std::vector<methods::ExecutionIssue> issues;
};

// Result of the application-facing facade, which performs applicability and automatic selection.
struct ApplicationCalculationResult {
    std::optional<charges::ChargeSet> charges;
    methods::ApplicabilityResult applicability;
    std::optional<ExecutionPolicy> execution_policy;
    std::vector<methods::ExecutionIssue> execution_issues;
    std::optional<methods::MethodOptions> effective_method_options;
    CalculationMetrics metrics;
    bool cancelled = false;

    [[nodiscard]] auto calculated() const noexcept -> bool {
        return charges.has_value();
    }
};

// Result of application-facing applicability and execution-plan assessment. The contained
// candidates reference methods from the registry and the owned parameter sets below.
struct ApplicationAssessmentResult {
    std::vector<parameters::ParameterSet> parameter_sets;
    methods::ApplicabilityResult applicability;
    const methods::ApplicableMethod* selected = nullptr;
    std::optional<ExecutionPolicy> execution_policy;
    std::vector<methods::ExecutionIssue> execution_issues;

    [[nodiscard]] auto executable() const noexcept -> bool {
        return execution_policy.has_value();
    }
};

// Owns facade inputs so language bindings and adapters do not need to manage native method
// pointers, parameter spans, or prepared-feature lifetimes. This request performs applicability
// using classification_options, chooses the highest-ranked applicable candidate, then executes it.
struct ApplicationCalculationRequest {
    core::MoleculeCollection molecules;
    std::vector<parameters::ParameterSet> parameter_sets;
    std::optional<std::string> method_id;
    std::optional<std::string> parameter_set_id;
    std::unordered_map<std::string, methods::MethodOptions> method_options;
    parameters::ClassificationOptions classification_options{};
    ExecutionSelection execution_selection{};
    ResourcePolicy resource_policy{};
    // Optional, non-owning observer for progress and cancellation. Nullptr means no observation.
    const CalculationObserver* observer = nullptr;
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

// Executes an already applicable and selected candidate. Its stored parameter classifications are
// used directly; classification policy is not reconsidered during execution.
[[nodiscard]] auto calculate(const CalculationRequest& request) -> CalculationResult;

// Convenience facade assessment: resolves registered methods, applies classification policy, and
// selects an execution plan without running a charge calculation. Explicit unavailable method or
// parameter-set IDs are errors; an unavailable explicit execution plan is reported as no plan.
[[nodiscard]] auto assess(const ApplicationCalculationRequest& request)
    -> ApplicationAssessmentResult;

// Convenience facade: calculates an owned molecule collection using registered methods and caller-
// owned parameter data. Omitted IDs enable deterministic automatic selection. Specified IDs
// restrict selection to the exact method or parameter set; unavailable or inapplicable selections
// are errors.
[[nodiscard]] auto calculate(const ApplicationCalculationRequest& request)
    -> ApplicationCalculationResult;

} // namespace chargefw::calculation
