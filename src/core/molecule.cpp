#include "chargefw/core/molecule.h"

#include "chargefw/core/atom.h"
#include "chargefw/core/bond.h"

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string_view>
#include <algorithm>

namespace chargefw::core {

namespace {
auto validate_conformers(const std::vector<Atom>& atoms, const std::vector<Conformer>& conformers)
    -> void {
    const auto atom_count = atoms.size();

    for (const auto& conformer : conformers) {
        if (conformer.size() != atom_count) {
            throw std::invalid_argument{"conformer position count must match molecule atom count"};
        }
    }
}

auto validate_bond_atom_indices(
    const std::vector<Atom>& atoms,
    const std::vector<Bond>& bonds
) -> void
{
    const auto atom_count = atoms.size();

    for (const auto& bond : bonds) {
        if (
            bond.first_atom_index() >= atom_count ||
            bond.second_atom_index() >= atom_count
        ) {
            throw std::invalid_argument{"bond refers to an atom index outside the molecule"};
        }
    }
}

auto same_unordered_bond(const Bond& first, const Bond& second) noexcept -> bool
{
    return (
        first.first_atom_index() == second.first_atom_index() &&
        first.second_atom_index() == second.second_atom_index()
    ) || (
        first.first_atom_index() == second.second_atom_index() &&
        first.second_atom_index() == second.first_atom_index()
    );
}

auto validate_no_duplicate_bonds(const std::vector<Bond>& bonds) -> void
{
    for (auto first = bonds.begin(); first != bonds.end(); ++first) {
        const auto duplicate = std::find_if(
            std::next(first),
            bonds.end(),
            [first](const Bond& second) -> bool {
                return same_unordered_bond(*first, second);
            }
        );

        if (duplicate != bonds.end()) {
            throw std::invalid_argument{"molecule contains duplicate bonds"};
        }
    }
}


} // namespace

Molecule::Molecule(std::vector<Atom> atoms, std::vector<Bond> bonds,
                   std::vector<Conformer> conformers, std::string name)
    : atoms_{std::move(atoms)}, bonds_{std::move(bonds)}, conformers_{std::move(conformers)},
      name_{std::move(name)} {
    validate_conformers(atoms_, conformers_);
    validate_bond_atom_indices(atoms_, bonds_);
    validate_no_duplicate_bonds(bonds_);
}

auto Molecule::name() const noexcept -> std::string_view {
    return name_;
}

auto Molecule::atoms() const noexcept -> std::span<const Atom> {
    return atoms_;
}

auto Molecule::bonds() const noexcept -> std::span<const Bond> {
    return bonds_;
}

auto Molecule::conformers() const noexcept -> std::span<const Conformer> {
    return conformers_;
}

auto Molecule::atom_count() const noexcept -> std::size_t {
    return atoms_.size();
}

auto Molecule::bond_count() const noexcept -> std::size_t {
    return bonds_.size();
}

auto Molecule::conformer_count() const noexcept -> std::size_t {
    return conformers_.size();
}

auto Molecule::atom(const std::size_t index) const -> const Atom& {
    return atoms_.at(index);
}

auto Molecule::bond(const std::size_t index) const -> const Bond& {
    return bonds_.at(index);
}

auto Molecule::conformer(const std::size_t index) const -> const Conformer& {
    return conformers_.at(index);
}

auto Molecule::has_coordinates() const noexcept -> bool {
    return !conformers_.empty();
}

} // namespace chargefw::core
