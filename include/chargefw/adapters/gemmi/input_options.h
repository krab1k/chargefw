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
};

} // namespace chargefw::adapters::gemmi
