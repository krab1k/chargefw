#include "support/test_molecules.h"

#include <chargefw/core/molecule_collection.h>

#include <stdexcept>
#include <string_view>
#include <vector>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;

TEST_CASE("molecule collection preserves ordered molecules", "[core][molecule-collection]") {
    const core::MoleculeCollection collection{
        {chargefw::test::make_water(), chargefw::test::make_formally_charged_pair()}, "examples"};

    CHECK(collection.name() == std::string_view{"examples"});
    CHECK(collection.size() == 2);
    CHECK_FALSE(collection.empty());
    CHECK(collection.molecules().size() == 2);
    CHECK(collection[0].name() == std::string_view{"water"});
    CHECK(collection.at(1).name() == std::string_view{"charged-pair"});
}

TEST_CASE("molecule collection supports no molecules", "[core][molecule-collection]") {
    const core::MoleculeCollection collection{std::vector<core::Molecule>{}, "empty"};

    CHECK(collection.empty());
}

TEST_CASE("molecule collection bounds-checked access rejects an invalid index",
          "[core][molecule-collection]") {
    const core::MoleculeCollection collection{
        {chargefw::test::make_water(), chargefw::test::make_formally_charged_pair()}, "examples"};

    CHECK_THROWS_AS(collection.at(2), std::out_of_range);
}
