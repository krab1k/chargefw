#pragma once

#include <cstdint>

namespace chargefw::adapters::gemmi {

enum class RecordSelection : std::uint8_t {
    all,
    polymers_and_ligands,
    polymers,
};

enum class BondStrategy : std::uint8_t {
    none,
    explicit_bonds,
    templates,
    hybrid,
};

struct InputOptions {
    RecordSelection selection = RecordSelection::all;
    BondStrategy bond_strategy = BondStrategy::none;
};

} // namespace chargefw::adapters::gemmi
