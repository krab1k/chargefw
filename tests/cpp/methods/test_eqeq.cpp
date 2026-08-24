#include "support/test_calculation.h"
#include "support/test_molecules.h"

#include <cmath>
#include <snitch/snitch.hpp>

TEST_CASE("EQeq has a stable water regression and responds to geometry", "[methods][eqeq]") {
    const auto charge_set =
        chargefw::test::calculate_method(chargefw::test::make_two_conformer_water(), "eqeq");
    const auto& charges = charge_set.assignment(0).charges;

    CHECK(std::abs(charges[0] - (-0.36761156)) < 1.0e-8);
    CHECK(std::abs(charges[1] - (0.18380578)) < 1.0e-8);
    CHECK(std::abs(charges[2] - (0.18380578)) < 1.0e-8);
    CHECK(std::abs(charges[0] - charge_set.assignment(1).charges[0]) > 1.0e-8);
}
