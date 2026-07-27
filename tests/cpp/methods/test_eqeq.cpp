#include "support/test_calculation.h"
#include "support/test_molecules.h"

#include <cassert>
#include <cmath>

namespace {

auto assert_close(const double actual, const double expected) -> void {
    assert(std::abs(actual - expected) < 1.0e-8);
}

} // namespace

auto main() -> int {
    const auto charge_set =
        chargefw::test::calculate_method(chargefw::test::make_two_conformer_water(), "eqeq");
    const auto& charges = charge_set.assignment(0).charges;

    chargefw::test::assert_calculation_provenance(charge_set, "eqeq", std::nullopt);
    chargefw::test::assert_conformer_dependent(charge_set, 2);

    assert(charges.size() == 3);
    assert_close(charges[0], -0.36751024);
    assert_close(charges[1], 0.18377329);
    assert_close(charges[2], 0.18373695);
    assert_close(charges.total(), 0.0);
    assert(std::abs(charges[0] - charge_set.assignment(1).charges[0]) > 1.0e-8);

    return 0;
}
