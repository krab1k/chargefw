#include "support/test_calculation.h"
#include "support/test_molecules.h"

#include <cmath>
#include <snitch/snitch.hpp>

TEST_CASE("MGC has a stable water graph-charge regression", "[methods][mgc]") {
    const auto charge_set =
        chargefw::test::calculate_single_method(chargefw::test::make_water_graph(), "mgc");
    const auto& charges = charge_set.assignment(0).charges;

    CHECK(std::abs(charges[0] - (-0.24280469939476498)) < 1.0e-12);
    CHECK(std::abs(charges[1] - 0.12140234969738249) < 1.0e-12);
    CHECK(std::abs(charges[2] - 0.12140234969738249) < 1.0e-12);
}
