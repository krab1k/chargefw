#pragma once

#include <cstddef>
#include <limits>
#include <stdexcept>

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

} // namespace chargefw::calculation::detail
