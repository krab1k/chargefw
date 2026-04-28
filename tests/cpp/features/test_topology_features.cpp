#include "support/test_molecules.h"

#include <chargefw/features/topology_features.h>

#include <algorithm>
#include <cassert>
#include <span>
#include <stdexcept>
#include <vector>

namespace features = chargefw::features;

namespace {

auto to_vector(std::span<const std::size_t> values) -> std::vector<std::size_t> {
    return {values.begin(), values.end()};
}

} // namespace

auto main() -> int {
    const auto water = chargefw::test::make_water();
    const features::TopologyFeatures topology{water};

    assert(&topology.molecule() == &water);

    assert(topology.degree(0) == 2);
    assert(topology.degree(1) == 1);
    assert(topology.degree(2) == 1);

    auto oxygen_neighbors = to_vector(topology.neighbor_indices(0));
    std::ranges::sort(oxygen_neighbors);

    assert((oxygen_neighbors == std::vector<std::size_t>{1, 2}));

    const auto h1_neighbors = to_vector(topology.neighbor_indices(1));
    assert(h1_neighbors == std::vector<std::size_t>{0});

    const auto oh1_bond_index = topology.bond_index_between(0, 1);
    const auto oh2_bond_index = topology.bond_index_between(0, 2);
    const auto no_self_bond = topology.bond_index_between(0, 0);

    assert(oh1_bond_index.has_value());
    assert(oh2_bond_index.has_value());
    assert(!no_self_bond.has_value());

    assert(topology.are_bonded(0, 1));
    assert(topology.are_bonded(1, 0));
    assert(!topology.are_bonded(1, 2));

    const auto oxygen_bonds = to_vector(topology.incident_bond_indices(0));
    assert(oxygen_bonds.size() == 2);

    bool rejected_atom_index = false;

    try {
        [[maybe_unused]] const auto invalid_degree = topology.degree(3);
    } catch (const std::out_of_range&) {
        rejected_atom_index = true;
    }

    assert(rejected_atom_index);

    return 0;
}
