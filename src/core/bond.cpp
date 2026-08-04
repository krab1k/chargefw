#include <chargefw/core/bond.h>

#include <stdexcept>

namespace chargefw::core {
namespace {

auto validate_bond_indices(const std::size_t first_atom_index, const std::size_t second_atom_index)
    -> void {
    if (first_atom_index == second_atom_index) {
        throw std::invalid_argument{"bond endpoints must refer to different atoms"};
    }
}

} // namespace

auto bond_order_from_value(const int value) -> BondOrder {
    switch (value) {
    case 1:
        return BondOrder::SINGLE;
    case 2:
        return BondOrder::DOUBLE;
    case 3:
        return BondOrder::TRIPLE;
    default:
        throw std::invalid_argument{"unsupported bond order " + std::to_string(value)};
    }
}

Bond::Bond(const std::size_t first_atom_index, const std::size_t second_atom_index,
           const BondOrder order)
    : first_atom_index_{first_atom_index}, second_atom_index_{second_atom_index}, order_{order} {
    validate_bond_indices(first_atom_index_, second_atom_index_);
}

auto Bond::first_atom_index() const noexcept -> std::size_t {
    return first_atom_index_;
}

auto Bond::second_atom_index() const noexcept -> std::size_t {
    return second_atom_index_;
}

auto Bond::order() const noexcept -> BondOrder {
    return order_;
}

} // namespace chargefw::core
