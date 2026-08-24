#include "support/test_molecules.h"

#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>

#include <array>
#include <future>
#include <type_traits>
#include <vector>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;
namespace features = chargefw::features;

static_assert(!std::is_constructible_v<features::PreparedMolecule, core::Molecule&&>);
static_assert(
    !std::is_constructible_v<features::PreparedMoleculeCollection, core::MoleculeCollection&&>);

TEST_CASE("prepared molecules retain source identity and cached topology", "[features][prepared]") {
    const auto water = chargefw::test::make_water_graph();
    const features::PreparedMolecule prepared{water};

    CHECK(&prepared.molecule() == &water);
    CHECK(&prepared.topology().molecule() == &water);
    CHECK(prepared.topology().degree(0) == 2);
}

TEST_CASE("prepared molecule collection retains each source molecule", "[features][prepared]") {
    const auto collection = core::MoleculeCollection{
        std::vector{chargefw::test::make_water_graph(), chargefw::test::make_water_graph()},
        "waters"};
    const features::PreparedMoleculeCollection prepared{collection};

    REQUIRE(prepared.size() == collection.size());
    CHECK(&prepared[0].molecule() == &collection[0]);
    CHECK(&prepared[1].molecule() == &collection[1]);
    CHECK(prepared[0].topology().degree(0) == 2);
    CHECK(prepared[1].topology().degree(0) == 2);
}

TEST_CASE("prepared topology supports concurrent read-only access", "[features][prepared]") {
    const auto water = chargefw::test::make_water_graph();
    const features::PreparedMolecule prepared{water};

    const auto read_topology = [&prepared] {
        for (auto iteration = 0; iteration < 1000; ++iteration) {
            if (prepared.topology().degree(0) != 2 || !prepared.topology().are_bonded(0, 1) ||
                !prepared.topology().are_bonded(0, 2)) {
                return false;
            }
        }

        return true;
    };

    std::array futures{std::async(std::launch::async, read_topology),
                       std::async(std::launch::async, read_topology),
                       std::async(std::launch::async, read_topology),
                       std::async(std::launch::async, read_topology)};

    for (auto& future : futures) {
        CHECK(future.get());
    }
}
