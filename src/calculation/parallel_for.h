#pragma once

#include <chargefw/calculation/observer.h>

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <utility>

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>

namespace chargefw::calculation::detail {

// A zero thread count delegates the worker count to oneTBB. The callable must write only to the
// result associated with its index; execution order is intentionally unspecified.
template <typename Function>
auto parallel_for_indexed(const std::size_t count, const std::size_t max_threads,
                          Function&& function) -> void {
    if (count == 0) {
        return;
    }

    if (max_threads == 1) {
        for (std::size_t index = 0; index < count; ++index) {
            function(index);
        }
        return;
    }

    if (max_threads > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument{"maximum thread count exceeds oneTBB's supported range"};
    }

    oneapi::tbb::task_arena arena{static_cast<int>(max_threads)};
    arena.execute([&] {
        oneapi::tbb::parallel_for(
            oneapi::tbb::blocked_range<std::size_t>{0, count},
            [&function](const oneapi::tbb::blocked_range<std::size_t>& range) {
                for (auto index = range.begin(); index != range.end(); ++index) {
                    function(index);
                }
            });
    });
}

// Context carried through a parallel loop for progress emission. observer, mode, and method_id are
// shared; target_index/molecule_index/conformer_index identify the enclosing target.
struct ProgressContext {
    const CalculationObserver* observer = nullptr;
    ExecutionMode mode{};
    std::string_view method_id;
    std::size_t target_index{};
    std::size_t target_count{};
    std::size_t molecule_index{};
    std::optional<std::size_t> conformer_index{};
    std::chrono::steady_clock::time_point computation_start{};
};

// Throws CalculationCancelled if the observer requests termination. Called by executors at
// cancellation check points (start of each target iteration).
inline auto check_cancellation(const CalculationObserver* observer) -> void {
    if (observer != nullptr && observer->cancelled()) {
        throw CalculationCancelled{};
    }
}

// Emits a single target-tier progress event. Called by executors at the start and end of each
// target iteration, where the per-target molecule_index and conformer_index are known.
inline auto emit_target_event(const CalculationObserver* observer, const CalculationPhase phase,
                              const ProgressContext& ctx) -> void {
    if (observer == nullptr) {
        return;
    }
    observer->on_progress(CalculationProgress{
        .phase = phase,
        .mode = ctx.mode,
        .method_id = ctx.method_id,
        .target_index = ctx.target_index,
        .target_count = ctx.target_count,
        .molecule_index = ctx.molecule_index,
        .conformer_index = ctx.conformer_index,
        .elapsed_seconds =
            std::chrono::duration<double>{std::chrono::steady_clock::now() - ctx.computation_start}
                .count(),
    });
}

// Minimum interval between fragment-tier progress events, in nanoseconds.
inline constexpr std::int64_t fragment_throttle_ns = 100'000'000; // 100 ms

// Fragment-tier loop with throttled progress emission and cooperative cancellation. When observer
// is null, delegates directly to parallel_for_indexed with zero overhead. Otherwise checks
// cancellation at the start of each iteration and emits a throttled fragment_finished event
// carrying the completed count, then a final event after the loop so the observer sees 100%
// completion.
//
// Target-tier events are NOT emitted here — executors emit them via check_cancellation() and
// emit_target_event() so per-target molecule/conformer identity is correctly populated.
template <typename Function>
auto progress_for_indexed(const std::size_t count, const std::size_t max_threads,
                          const ProgressContext& ctx, Function&& function) -> void {
    if (ctx.observer == nullptr || count == 0) {
        parallel_for_indexed(count, max_threads, std::forward<Function>(function));
        return;
    }

    const auto* observer = ctx.observer;

    const auto elapsed = [&ctx]() -> double {
        return std::chrono::duration<double>{std::chrono::steady_clock::now() -
                                             ctx.computation_start}
            .count();
    };

    std::atomic<std::int64_t> last_emit_ns{0};
    std::atomic<std::size_t> completed{0};

    auto wrapped = [&](const std::size_t index) {
        if (observer->cancelled()) {
            throw CalculationCancelled{};
        }
        function(index);
        const auto done = completed.fetch_add(1, std::memory_order_relaxed) + 1;

        const auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::steady_clock::now().time_since_epoch())
                                .count();
        auto last = last_emit_ns.load(std::memory_order_relaxed);
        if (now_ns - last >= fragment_throttle_ns) {
            if (last_emit_ns.compare_exchange_strong(last, now_ns, std::memory_order_relaxed)) {
                observer->on_progress(CalculationProgress{
                    .phase = CalculationPhase::fragment_finished,
                    .mode = ctx.mode,
                    .method_id = ctx.method_id,
                    .target_index = ctx.target_index,
                    .target_count = ctx.target_count,
                    .fragment_index = done,
                    .fragment_count = count,
                    .molecule_index = ctx.molecule_index,
                    .conformer_index = ctx.conformer_index,
                    .elapsed_seconds = elapsed(),
                });
            }
        }
    };
    parallel_for_indexed(count, max_threads, wrapped);

    // Final event so the observer sees 100% completion.
    observer->on_progress(CalculationProgress{
        .phase = CalculationPhase::fragment_finished,
        .mode = ctx.mode,
        .method_id = ctx.method_id,
        .target_index = ctx.target_index,
        .target_count = ctx.target_count,
        .fragment_index = count,
        .fragment_count = count,
        .molecule_index = ctx.molecule_index,
        .conformer_index = ctx.conformer_index,
        .elapsed_seconds = elapsed(),
    });
}

} // namespace chargefw::calculation::detail
