#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/core/periodic_table.h>

#include <cstddef>
#include <string>

namespace chargefw::core::detail {

[[nodiscard]] inline auto bond_order_description(const BondOrder order) -> std::string_view {
    switch (order) {
    case BondOrder::SINGLE:
        return "single";
    case BondOrder::DOUBLE:
        return "double";
    case BondOrder::TRIPLE:
        return "triple";
    }
    return "unknown";
}

[[nodiscard]] inline auto atom_description(const Molecule& molecule, const std::size_t atom_index)
    -> std::string {
    const auto& atom = molecule.atom(atom_index);
    auto description = "atom " + std::to_string(atom_index + 1) + " (";
    if (!atom.name().empty()) {
        description += "source name '" + std::string{atom.name()} + "', ";
    }
    description += std::string{element_symbol(atom.atomic_number())} + ", formal charge " +
                   (atom.formal_charge() > 0 ? "+" : "") + std::to_string(atom.formal_charge()) +
                   ")";
    return description;
}

[[nodiscard]] inline auto bond_description(const Molecule& molecule, const std::size_t bond_index)
    -> std::string {
    const auto& bond = molecule.bond(bond_index);
    return "bond " + std::to_string(bond_index + 1) + " between " +
           atom_description(molecule, bond.first_atom_index()) + " and " +
           atom_description(molecule, bond.second_atom_index()) + ", " +
           std::string{bond_order_description(bond.order())} + " bond";
}

[[nodiscard]] inline auto molecule_description(const Molecule& molecule,
                                               const std::size_t molecule_index) -> std::string {
    auto description = "molecule " + std::to_string(molecule_index + 1);
    if (!molecule.name().empty()) {
        description += " ('" + std::string{molecule.name()} + "')";
    }
    return description;
}

} // namespace chargefw::core::detail
