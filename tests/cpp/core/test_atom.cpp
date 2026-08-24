#include <chargefw/core/atom.h>

#include <stdexcept>
#include <string_view>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;

TEST_CASE("atom stores supplied properties", "[core][atom]") {
    const core::Atom carbon{6, 0, "C1"};

    CHECK(carbon.atomic_number() == 6);
    CHECK(carbon.formal_charge() == 0);
    CHECK(carbon.name() == std::string_view{"C1"});

    const core::Atom ammonium_nitrogen{7, 1, "N"};

    CHECK(ammonium_nitrogen.atomic_number() == 7);
    CHECK(ammonium_nitrogen.formal_charge() == 1);
    CHECK(ammonium_nitrogen.name() == std::string_view{"N"});

    const core::Atom unnamed_oxygen{8};

    CHECK(unnamed_oxygen.atomic_number() == 8);
    CHECK(unnamed_oxygen.formal_charge() == 0);
    CHECK(unnamed_oxygen.name().empty());
}

TEST_CASE("atom rejects invalid atomic numbers", "[core][atom]") {
    CHECK_THROWS_AS(core::Atom{0}, std::invalid_argument);
    CHECK_THROWS_AS(core::Atom{119}, std::invalid_argument);
}
