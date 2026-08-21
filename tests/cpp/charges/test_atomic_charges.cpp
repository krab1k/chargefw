#include "support/test_assertions.h"

#include <chargefw/charges/atomic_charges.h>

#include <cassert>
#include <limits>
#include <stdexcept>
#include <vector>

namespace charges = chargefw::charges;

auto main() -> int {
    const charges::AtomicCharges values{{0.2, -0.1, -0.1}};

    assert(values.size() == 3);
    assert(!values.empty());
    assert(values.values().size() == 3);
    chargefw::test::assert_close(values.total(), 0.0, 1e-12);
    chargefw::test::assert_close(values[0], 0.2, 1e-12);
    chargefw::test::assert_close(values.at(2), -0.1, 1e-12);

    const charges::AtomicCharges empty_values{std::vector<double>{}};
    assert(empty_values.empty());
    chargefw::test::assert_close(empty_values.total(), 0.0, 1e-12);

    bool rejected_bad_index = false;

    try {
        [[maybe_unused]] const auto value = values.at(3);
    } catch (const std::out_of_range&) {
        rejected_bad_index = true;
    }

    assert(rejected_bad_index);

    bool rejected_non_finite = false;

    try {
        [[maybe_unused]] const charges::AtomicCharges invalid{
            {0.0, std::numeric_limits<double>::infinity()}};
    } catch (const std::invalid_argument&) {
        rejected_non_finite = true;
    }

    assert(rejected_non_finite);

    return 0;
}
