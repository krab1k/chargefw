#pragma once

#include <chargefw/calculation/assessment.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>

#include <optional>
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
struct ApplicationExecutionResult {
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

// Executes an already applicable and selected candidate. Its stored parameter classifications are
// used directly; classification policy is not reconsidered during execution.
[[nodiscard]] auto calculate(const CalculationRequest& request) -> CalculationResult;

// Executes a plan returned by assess() without repeating preparation or applicability assessment.
// The assessment is consumed because the result owns its classifications and parameter data.
[[nodiscard]] auto calculate(ApplicationAssessmentResult assessment, std::size_t max_threads = 0,
                             const CalculationObserver* observer = nullptr)
    -> ApplicationExecutionResult;

} // namespace chargefw::calculation
