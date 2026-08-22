#pragma once

#include <chargefw/calculation/assessment.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
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

// Owns facade inputs so language bindings and adapters do not need to manage native method
// pointers, parameter spans, or prepared-feature lifetimes. This request performs applicability,
// selection, and execution through the one-shot calculate() convenience facade.
struct ApplicationCalculationRequest {
    core::MoleculeCollection molecules;
    std::vector<parameters::ParameterSet> parameter_sets;
    std::optional<std::string> method_id;
    std::optional<std::string> parameter_set_id;
    std::unordered_map<std::string, methods::MethodOptions> method_options;
    parameters::ClassificationOptions classification_options{};
    ExecutionSelection execution_selection{};
    ResourcePolicy resource_policy{};
    // Optional, non-owning observer for the one-shot calculation convenience facade.
    const CalculationObserver* observer = nullptr;
};

// Executes an already applicable and selected candidate. Its stored parameter classifications are
// used directly; classification policy is not reconsidered during execution.
[[nodiscard]] auto calculate(const CalculationRequest& request) -> CalculationResult;

// Executes a plan returned by assess() without repeating preparation or applicability assessment.
// The assessment is consumed because the result owns its classifications and parameter data.
[[nodiscard]] auto calculate(ApplicationAssessmentResult assessment, std::size_t max_threads = 0,
                             const CalculationObserver* observer = nullptr)
    -> ApplicationCalculationResult;

// Convenience facade: calculates an owned molecule collection using registered methods and caller-
// owned parameter data. Omitted IDs enable deterministic automatic selection. Specified IDs
// restrict selection to the exact method or parameter set; unavailable or inapplicable selections
// are errors.
[[nodiscard]] auto calculate(const ApplicationCalculationRequest& request)
    -> ApplicationCalculationResult;

} // namespace chargefw::calculation
