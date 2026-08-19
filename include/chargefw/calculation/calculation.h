#pragma once

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/parameters/classification/classification_result.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <optional>
#include <string>
#include <vector>

namespace chargefw::calculation {

struct CalculationRequest {
    // Non-owning input; the referenced collection and selected candidate must outlive calculate().
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
    const features::PreparedMoleculeCollection& molecules;
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
    const methods::ApplicableMethod& selected;
};

struct CalculationResult {
    charges::ChargeSet charges;
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

    [[nodiscard]] auto calculated() const noexcept -> bool {
        return charges.has_value();
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
    parameters::ClassificationOptions classification_options{};
    ExecutionSelection execution_selection{};
    ResourcePolicy resource_policy{};
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

// Convenience facade: calculates an owned molecule collection using registered methods and caller-
// owned parameter data. Omitted IDs enable deterministic automatic selection. Specified IDs
// restrict selection to the exact method or parameter set; unavailable or inapplicable selections
// are errors.
[[nodiscard]] auto calculate(const ApplicationCalculationRequest& request)
    -> ApplicationCalculationResult;

} // namespace chargefw::calculation
