#pragma once

namespace chargefw::adapters::gemmi {

enum class RecordSelection {
    all,
    polymers_and_ligands,
    polymers,
};

enum class BondStrategy {
    none,
    explicit_bonds,
    templates,
    hybrid,
};

} // namespace chargefw::adapters::gemmi
