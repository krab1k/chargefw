#include "support/test_calculation.h"
#include "support/test_molecules.h"

#include <cmath>
#include <snitch/snitch.hpp>

TEST_CASE("EQeq responds to geometry", "[methods][eqeq]") {
    const auto charge_set =
        chargefw::test::calculate_method(chargefw::test::make_two_conformer_water(), "eqeq");

    CHECK(std::abs(charge_set.assignment(0).charges[0] - charge_set.assignment(1).charges[0]) >
          1.0e-8);
}
