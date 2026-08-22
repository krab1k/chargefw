#pragma once

#include <chargefw/calculation/execution_policy.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string_view>

namespace chargefw::calculation {

// Coarse calculation phases reported by the facade and executors. Target-level events are always
// emitted (rare); fragment-level events are throttled by progress_for_indexed to avoid starving the
// solver on large inputs.
enum class CalculationPhase : std::uint8_t {
    // Facade-level computation phases, emitted once each.
    computation_started,
    computation_finished,
    // Target tier: one molecule+conformer pair.
    target_started,
    target_finished,
    // Fragment tier: one center-atom (cutoff) or pivot (cover) fragment.
    fragment_started,
    fragment_finished,
};

// Structured progress snapshot. Fields are only populated for the phase that produced the event.
// The non-owning method_id view is valid only during the on_progress callback; callers that need to
// retain it must copy into owned storage.
struct CalculationProgress {
    CalculationPhase phase{};

    // Effective execution mode of the running calculation.
    ExecutionMode mode{};

    // Selected method ID; valid from computation_started onward.
    std::string_view method_id{};

    // Target tier: current and total target indices across all molecules/conformers.
    std::size_t target_index{};
    std::size_t target_count{};

    // Fragment tier: current and total fragment indices within the current target.
    std::size_t fragment_index{};
    std::size_t fragment_count{};

    // Molecule and conformer identifying the current target.
    std::size_t molecule_index{};
    std::optional<std::size_t> conformer_index{};

    // Elapsed computation seconds at the time of the event.
    double elapsed_seconds{};
};

// Cooperative, per-request observer of calculation progress. The observer is non-owning: callers
// retain the object and pass a reference through calculation requests. Executor callbacks can run
// from oneTBB worker threads; implementations must be thread-safe.
//
// The observer is purely observational: it must not mutate method options, parameters, execution
// policy, geometry, or selection. Implementations must not throw; exceptions would propagate as
// calculation failures. Returning true from cancelled() requests early termination; the application
// facade converts this into a cancelled result with no partial charges.
class CalculationObserver {
  public:
    CalculationObserver() = default;
    CalculationObserver(const CalculationObserver&) = delete;
    auto operator=(const CalculationObserver&) -> CalculationObserver& = delete;
    CalculationObserver(CalculationObserver&&) = delete;
    auto operator=(CalculationObserver&&) -> CalculationObserver& = delete;
    virtual ~CalculationObserver() = default;

    virtual void on_progress(const CalculationProgress& /*progress*/) const {}

    // When true, the running calculation is asked to stop as soon as control returns to a
    // cancellation check point. Checked at the fragment tier inside progress_for_indexed.
    [[nodiscard]] virtual auto cancelled() const noexcept -> bool {
        return false;
    }
};

// Stateless observer used when callers do not need progress or cancellation handling.
[[nodiscard]] inline auto default_calculation_observer() -> const CalculationObserver& {
    static const auto observer = CalculationObserver{};
    return observer;
}

// Thrown when a cooperative cancellation request is observed mid-calculation. The facade catches
// it and converts the result to a clean cancellation rather than a partial charge set.
class CalculationCancelled : public std::runtime_error {
  public:
    CalculationCancelled() : std::runtime_error{"calculation was cancelled"} {}
};

} // namespace chargefw::calculation
