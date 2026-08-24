#include "support/test_molecules.h"

#include <chargefw/features/topology_features.h>

#include <algorithm>
#include <span>
#include <stdexcept>
#include <vector>

#include <snitch/snitch.hpp>

namespace features = chargefw::features;

namespace {

auto to_vector(std::span<const std::size_t> values) -> std::vector<std::size_t> {
    return {values.begin(), values.end()};
}

} // namespace

TEST_CASE("topology features expose molecule reference and graph degree", "[features][topology]") {
    const auto water = chargefw::test::make_water();
    const features::TopologyFeatures topology{water};

    CHECK(&topology.molecule() == &water);
    CHECK(topology.degree(0) == 2);
    CHECK(topology.degree(1) == 1);
    CHECK(topology.degree(2) == 1);
}

TEST_CASE("topology features expose sorted neighbor indices", "[features][topology]") {
    const auto water = chargefw::test::make_water();
    const features::TopologyFeatures topology{water};

    auto oxygen_neighbors = to_vector(topology.neighbor_indices(0));
    std::ranges::sort(oxygen_neighbors);
    CHECK(oxygen_neighbors == std::vector<std::size_t>{1, 2});

    const auto h1_neighbors = to_vector(topology.neighbor_indices(1));
    CHECK(h1_neighbors == std::vector<std::size_t>{0});
}

TEST_CASE("topology features resolve bond indices and adjacency", "[features][topology]") {
    const auto water = chargefw::test::make_water();
    const features::TopologyFeatures topology{water};

    const auto oh1_bond_index = topology.bond_index_between(0, 1);
    const auto oh2_bond_index = topology.bond_index_between(0, 2);
    const auto no_self_bond = topology.bond_index_between(0, 0);

    CHECK(oh1_bond_index.has_value());
    CHECK(oh2_bond_index.has_value());
    CHECK_FALSE(no_self_bond.has_value());

    CHECK(topology.are_bonded(0, 1));
    CHECK(topology.are_bonded(1, 0));
    CHECK_FALSE(topology.are_bonded(1, 2));

    const auto oxygen_bonds = to_vector(topology.incident_bond_indices(0));
    CHECK(oxygen_bonds.size() == 2);
}

TEST_CASE("topology features reject invalid atom indices", "[features][topology]") {
    const auto water = chargefw::test::make_water();
    const features::TopologyFeatures topology{water};

    CHECK_THROWS_AS(topology.degree(3), std::out_of_range);
}
