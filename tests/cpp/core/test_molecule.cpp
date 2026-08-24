#include "support/test_molecules.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>

#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;

TEST_CASE("molecule exposes source-ordered graph and conformer data", "[core][molecule]") {
    const auto water = chargefw::test::make_water();

    CHECK(water.name() == std::string_view{"water"});
    REQUIRE(water.atom_count() == 3);
    REQUIRE(water.bond_count() == 2);
    REQUIRE(water.conformer_count() == 1);
    CHECK(water.has_coordinates());
    CHECK(water.atoms().size() == 3);
    CHECK(water.bonds().size() == 2);
    CHECK(water.conformers().size() == 1);
    CHECK(water.atom(0).atomic_number() == 8);
    CHECK(water.bond(0).order() == core::BondOrder::SINGLE);
    CHECK(water.conformer(0).size() == 3);
}

TEST_CASE("molecule reports absent coordinates", "[core][molecule]") {
    const auto charged_pair = chargefw::test::make_formally_charged_pair();

    CHECK_FALSE(charged_pair.has_coordinates());
}

TEST_CASE("molecule bounds-checked atom access rejects an invalid index", "[core][molecule]") {
    const auto water = chargefw::test::make_water();

    CHECK_THROWS_AS(water.atom(3), std::out_of_range);
}

TEST_CASE("molecule rejects bonds with invalid atom indices", "[core][molecule]") {
    const auto make_invalid_molecule = [] {
        std::vector atoms{core::Atom{6, 0, "C"}};
        std::vector bonds{core::Bond{0, 1, core::BondOrder::SINGLE}};
        return core::Molecule{std::move(atoms), std::move(bonds)};
    };

    CHECK_THROWS_AS(make_invalid_molecule(), std::invalid_argument);
}

TEST_CASE("molecule rejects duplicate bonds", "[core][molecule]") {
    const auto make_invalid_molecule = [] {
        std::vector atoms{core::Atom{6, 0, "C"}, core::Atom{1, 0, "H"}};
        std::vector bonds{core::Bond{0, 1, core::BondOrder::SINGLE},
                          core::Bond{1, 0, core::BondOrder::SINGLE}};
        return core::Molecule{std::move(atoms), std::move(bonds)};
    };

    CHECK_THROWS_AS(make_invalid_molecule(), std::invalid_argument);
}

TEST_CASE("molecule rejects conformers with a mismatched atom count", "[core][molecule]") {
    const auto make_invalid_molecule = [] {
        std::vector atoms{core::Atom{8, 0, "O"}, core::Atom{1, 0, "H1"}};
        std::vector positions{core::Position{.x = 0.0, .y = 0.0, .z = 0.0}};
        std::vector conformers{core::Conformer{std::move(positions)}};
        return core::Molecule{std::move(atoms), {}, std::move(conformers)};
    };

    CHECK_THROWS_AS(make_invalid_molecule(), std::invalid_argument);
}
