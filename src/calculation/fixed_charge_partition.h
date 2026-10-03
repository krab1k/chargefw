#pragma once

#include <chargefw/calculation/assessment.h>

#include <cstddef>
#include <string>
#include <vector>

namespace chargefw::calculation::detail {

struct FixedChargePartitionTarget {
    // Active indices map to atoms/bonds in the original molecule at this target index.
    std::vector<std::size_t> active_atom_indices;
    std::vector<std::size_t> active_bond_indices;
    std::vector<FixedAtomCharge> sources;
    // Affected targets use original conformer order and sorted source atom order.
    std::vector<std::vector<core::Position>> source_positions;
    double original_charge = 0.0;
    double active_charge = 0.0;
};

struct FixedChargePartition {
    core::MoleculeCollection active_molecules;
    std::vector<FixedChargePartitionTarget> targets;
    std::string charge_provenance;
};

auto validate_fixed_charge_embedding(const core::MoleculeCollection& molecules,
                                     const FixedChargeEmbedding& embedding) -> void;

[[nodiscard]] auto make_fixed_charge_partition(const core::MoleculeCollection& molecules,
                                               const FixedChargeEmbedding& embedding)
    -> FixedChargePartition;

} // namespace chargefw::calculation::detail
