#include <chargefw/features/spatial_fragment.h>

#include <chargefw/features/conformer_features.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::features {
namespace {

auto validate_radius(const double radius) -> void {
    if (!std::isfinite(radius) || radius <= 0.0) {
        throw std::invalid_argument{"fragment radius must be finite and strictly positive"};
    }
}

} // namespace

SpatialFragment::SpatialFragment(core::Molecule molecule,
                                 std::vector<std::size_t> local_to_source_atom_indices,
                                 std::vector<std::size_t> source_to_local_atom_indices,
                                 std::vector<std::size_t> local_to_source_bond_indices,
                                 const std::size_t center_local_atom_index,
                                 const std::size_t source_conformer_index)
    : molecule_{std::move(molecule)},
      local_to_source_atom_indices_{std::move(local_to_source_atom_indices)},
      source_to_local_atom_indices_{std::move(source_to_local_atom_indices)},
      local_to_source_bond_indices_{std::move(local_to_source_bond_indices)},
      center_local_atom_index_{center_local_atom_index},
      source_conformer_index_{source_conformer_index} {}

auto SpatialFragment::molecule() const noexcept -> const core::Molecule& {
    return molecule_;
}

auto SpatialFragment::source_conformer_index() const noexcept -> std::size_t {
    return source_conformer_index_;
}

auto SpatialFragment::center_local_atom_index() const noexcept -> std::size_t {
    return center_local_atom_index_;
}

auto SpatialFragment::local_to_source_atom_indices() const noexcept
    -> std::span<const std::size_t> {
    return local_to_source_atom_indices_;
}

auto SpatialFragment::source_to_local_atom_indices() const noexcept
    -> std::span<const std::size_t> {
    return source_to_local_atom_indices_;
}

auto SpatialFragment::local_to_source_bond_indices() const noexcept
    -> std::span<const std::size_t> {
    return local_to_source_bond_indices_;
}

auto build_spatial_fragment(const PreparedMolecule& source, const std::size_t conformer_index,
                            const std::size_t center_atom_index, const double radius)
    -> SpatialFragment {
    validate_radius(radius);

    const auto& molecule = source.molecule();
    const ConformerFeatures geometry{molecule, conformer_index};
    auto selected_source_atom_indices = geometry.neighbor_indices_within(center_atom_index, radius);
    selected_source_atom_indices.push_back(center_atom_index);
    std::ranges::sort(selected_source_atom_indices);

    std::vector<std::size_t> source_to_local_atom_indices(molecule.atom_count(), no_source_index);
    std::vector<core::Atom> atoms;
    std::vector<core::Position> positions;
    atoms.reserve(selected_source_atom_indices.size());
    positions.reserve(selected_source_atom_indices.size());

    auto center_local_atom_index = no_source_index;
    for (std::size_t local_index = 0; local_index < selected_source_atom_indices.size();
         ++local_index) {
        const auto source_index = selected_source_atom_indices[local_index];
        const auto& source_atom = molecule.atom(source_index);
        atoms.emplace_back(source_atom.atomic_number(), source_atom.formal_charge(),
                           std::string{source_atom.name()});
        positions.push_back(geometry.position(source_index));
        source_to_local_atom_indices[source_index] = local_index;

        if (source_index == center_atom_index) {
            center_local_atom_index = local_index;
        }
    }

    std::vector<core::Bond> bonds;
    std::vector<std::size_t> local_to_source_bond_indices;
    bonds.reserve(molecule.bond_count());
    local_to_source_bond_indices.reserve(molecule.bond_count());

    for (std::size_t source_bond_index = 0; source_bond_index < molecule.bond_count();
         ++source_bond_index) {
        const auto& source_bond = molecule.bond(source_bond_index);
        const auto first_local = source_to_local_atom_indices[source_bond.first_atom_index()];
        const auto second_local = source_to_local_atom_indices[source_bond.second_atom_index()];
        if (first_local == no_source_index || second_local == no_source_index) {
            continue;
        }

        bonds.emplace_back(first_local, second_local, source_bond.order());
        local_to_source_bond_indices.push_back(source_bond_index);
    }

    std::vector<core::Conformer> conformers;
    conformers.emplace_back(std::move(positions),
                            std::string{molecule.conformer(conformer_index).name()});
    return SpatialFragment{core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers),
                                          std::string{molecule.name()}},
                           std::move(selected_source_atom_indices),
                           std::move(source_to_local_atom_indices),
                           std::move(local_to_source_bond_indices),
                           center_local_atom_index,
                           conformer_index};
}

auto project_classification(const parameters::ParameterClassification& source,
                            const SpatialFragment& fragment)
    -> parameters::ParameterClassification {
    std::vector<std::size_t> atom_parameter_entry_indices;
    if (!source.atom().empty()) {
        atom_parameter_entry_indices.reserve(fragment.local_to_source_atom_indices().size());
        for (const auto source_atom_index : fragment.local_to_source_atom_indices()) {
            atom_parameter_entry_indices.push_back(source.atom().at(source_atom_index));
        }
    }

    std::vector<std::size_t> bond_parameter_entry_indices;
    if (!source.bond().empty()) {
        bond_parameter_entry_indices.reserve(fragment.local_to_source_bond_indices().size());
        for (const auto source_bond_index : fragment.local_to_source_bond_indices()) {
            bond_parameter_entry_indices.push_back(source.bond().at(source_bond_index));
        }
    }

    return parameters::ParameterClassification{
        parameters::AtomParameterClassification{std::move(atom_parameter_entry_indices)},
        parameters::BondParameterClassification{std::move(bond_parameter_entry_indices)}};
}

} // namespace chargefw::features
