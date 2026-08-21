#pragma once

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace chargefw::test {

template <typename Exception, typename Callable>
[[nodiscard]] auto throws(Callable&& callable) -> bool {
    try {
        std::forward<Callable>(callable)();
    } catch (const Exception&) {
        return true;
    }

    return false;
}

template <typename Callable>
[[nodiscard]] auto throws_invalid_argument(Callable&& callable) -> bool {
    return throws<std::invalid_argument>(std::forward<Callable>(callable));
}

inline auto assert_close(const double actual, const double expected, const double tolerance)
    -> void {
    assert(std::abs(actual - expected) < tolerance);
}

} // namespace chargefw::test
