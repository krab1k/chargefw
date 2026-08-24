#include <chargefw/charges/atomic_charges.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <snitch/snitch.hpp>

namespace charges = chargefw::charges;

TEST_CASE("atomic charges expose values and totals", "[charges][atomic-charges]") {
    const charges::AtomicCharges values{{0.2, -0.1, -0.1}};

    REQUIRE(values.size() == 3);
    CHECK_FALSE(values.empty());
    CHECK(values.values().size() == 3);
    CHECK(std::abs(values.total() - 0.0) < 1e-12);
    CHECK(std::abs(values[0] - 0.2) < 1e-12);
    CHECK(std::abs(values.at(2) - (-0.1)) < 1e-12);
}

TEST_CASE("atomic charges support an empty value list", "[charges][atomic-charges]") {
    const charges::AtomicCharges empty_values{std::vector<double>{}};

    CHECK(empty_values.empty());
    CHECK(std::abs(empty_values.total() - 0.0) < 1e-12);
}

TEST_CASE("atomic charges reject non-finite values and invalid indices",
          "[charges][atomic-charges]") {
    const charges::AtomicCharges values{{0.2, -0.1, -0.1}};

    CHECK_THROWS_AS(values.at(3), std::out_of_range);
    CHECK_THROWS_AS((charges::AtomicCharges{{0.0, std::numeric_limits<double>::infinity()}}),
                    std::invalid_argument);
}
