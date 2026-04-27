#include <chargefw/features/prepared_molecule.h>

#include <stdexcept>

namespace chargefw::features {
namespace {

auto validate_atom_index(const core::Molecule& molecule, const std::size_t atom_index) -> void {
    if (atom_index >= molecule.atom_count()) {
        throw std::out_of_range{"atom index is outside the molecule"};
    }
}

} // namespace

PreparedMolecule::PreparedMolecule(const core::Molecule& molecule)
    : molecule_{&molecule}, neighbor_indices_(molecule.atom_count()),
      incident_bond_indices_(molecule.atom_count()) {
    const auto bonds = molecule.bonds();

    for (std::size_t bond_index = 0; bond_index < bonds.size(); ++bond_index) {
        const auto& bond = bonds[bond_index];

        const auto first_atom_index = bond.first_atom_index();
        const auto second_atom_index = bond.second_atom_index();

        neighbor_indices_[first_atom_index].push_back(second_atom_index);
        incident_bond_indices_[first_atom_index].push_back(bond_index);

        neighbor_indices_[second_atom_index].push_back(first_atom_index);
        incident_bond_indices_[second_atom_index].push_back(bond_index);
    }
}

auto PreparedMolecule::molecule() const noexcept -> const core::Molecule& {
    return *molecule_;
}

auto PreparedMolecule::degree(const std::size_t atom_index) const -> std::size_t {
    validate_atom_index(molecule(), atom_index);
    return neighbor_indices_[atom_index].size();
}

auto PreparedMolecule::neighbor_indices(const std::size_t atom_index) const
    -> std::span<const std::size_t> {
    validate_atom_index(molecule(), atom_index);

    const auto& neighbors = neighbor_indices_[atom_index];
    return {neighbors.data(), neighbors.size()};
}

auto PreparedMolecule::incident_bond_indices(const std::size_t atom_index) const
    -> std::span<const std::size_t> {
    validate_atom_index(molecule(), atom_index);

    const auto& bonds = incident_bond_indices_[atom_index];
    return {bonds.data(), bonds.size()};
}

auto PreparedMolecule::bond_index_between(const std::size_t first_atom_index,
                                          const std::size_t second_atom_index) const
    -> std::optional<std::size_t> {
    validate_atom_index(molecule(), first_atom_index);
    validate_atom_index(molecule(), second_atom_index);

    if (first_atom_index == second_atom_index) {
        return std::nullopt;
    }

    const auto& neighbors = neighbor_indices_[first_atom_index];
    const auto& incident_bonds = incident_bond_indices_[first_atom_index];

    for (std::size_t index = 0; index < neighbors.size(); ++index) {
        if (neighbors[index] == second_atom_index) {
            return incident_bonds[index];
        }
    }

    return std::nullopt;
}

auto PreparedMolecule::are_bonded(const std::size_t first_atom_index,
                                  const std::size_t second_atom_index) const -> bool {
    return bond_index_between(first_atom_index, second_atom_index).has_value();
}

} // namespace chargefw::features