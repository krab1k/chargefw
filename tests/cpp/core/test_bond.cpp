#include <chargefw/core/bond.h>

#include <cassert>
#include <stdexcept>

namespace core = chargefw::core;

auto main() -> int {
    const core::Bond bond{0, 1, core::BondOrder::SINGLE};

    assert(bond.first_atom_index() == 0);
    assert(bond.second_atom_index() == 1);
    assert(bond.order() == core::BondOrder::SINGLE);

    const core::Bond double_bond{2, 3, core::BondOrder::DOUBLE};

    assert(double_bond.first_atom_index() == 2);
    assert(double_bond.second_atom_index() == 3);
    assert(double_bond.order() == core::BondOrder::DOUBLE);

    bool rejected_self_bond = false;

    try {
        [[maybe_unused]] const core::Bond invalid{0, 0, core::BondOrder::SINGLE};
    } catch (const std::invalid_argument&) {
        rejected_self_bond = true;
    }

    assert(rejected_self_bond);

    return 0;
}