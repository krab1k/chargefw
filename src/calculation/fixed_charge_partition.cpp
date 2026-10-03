#include "fixed_charge_partition.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::calculation::detail {
namespace {

struct ValidatedTarget {
    std::vector<FixedAtomCharge> sources;
    double original_charge = 0.0;
    double active_charge = 0.0;
};

[[nodiscard]] auto validate_and_group_sources(const core::MoleculeCollection& molecules,
                                              const FixedChargeEmbedding& embedding)
    -> std::vector<ValidatedTarget> {
    auto selected_sources = embedding.sources;
    for (const auto& source : embedding.sources) {
        if (source.molecule_index >= molecules.size()) {
            throw std::invalid_argument{
                "fixed charge source molecule index " + std::to_string(source.molecule_index) +
                " is out of range (molecule count " + std::to_string(molecules.size()) + ")"};
        }
        const auto& molecule = molecules[source.molecule_index];
        if (source.atom_index >= molecule.atom_count()) {
            throw std::invalid_argument{
                "fixed charge source at molecule " + std::to_string(source.molecule_index) +
                " has atom index " + std::to_string(source.atom_index) +
                " out of range (atom count " + std::to_string(molecule.atom_count()) + ")"};
        }
        if (!std::isfinite(source.charge)) {
            throw std::invalid_argument{
                "fixed charge source at molecule " + std::to_string(source.molecule_index) +
                ", atom " + std::to_string(source.atom_index) + " has non-finite charge"};
        }
    }

    std::ranges::sort(selected_sources, {}, [](const FixedAtomCharge& source) {
        return std::pair{source.molecule_index, source.atom_index};
    });
    for (std::size_t index = 1; index < selected_sources.size(); ++index) {
        if (selected_sources[index].molecule_index == selected_sources[index - 1].molecule_index &&
            selected_sources[index].atom_index == selected_sources[index - 1].atom_index) {
            throw std::invalid_argument{"duplicate fixed charge source at molecule " +
                                        std::to_string(selected_sources[index].molecule_index) +
                                        ", atom " +
                                        std::to_string(selected_sources[index].atom_index)};
        }
    }

    auto targets = std::vector<ValidatedTarget>(molecules.size());

    for (auto first = selected_sources.begin(); first != selected_sources.end();) {
        const auto molecule_index = first->molecule_index;
        const auto last = std::ranges::find_if(first, selected_sources.end(),
                                               [molecule_index](const FixedAtomCharge& selected) {
                                                   return selected.molecule_index != molecule_index;
                                               });
        const auto& molecule = molecules[molecule_index];
        auto& target = targets[molecule_index];
        target.sources.assign(first, last);

        auto selected_mask = std::vector<bool>(molecule.atom_count(), false);
        for (auto selected = first; selected != last; ++selected) {
            selected_mask[selected->atom_index] = true;
        }
        for (const auto& bond : molecule.bonds()) {
            if (selected_mask[bond.first_atom_index()]) {
                throw std::invalid_argument{
                    "fixed charge source at molecule " + std::to_string(molecule_index) +
                    ", atom " + std::to_string(bond.first_atom_index()) + " is bonded to atom " +
                    std::to_string(bond.second_atom_index())};
            }
            if (selected_mask[bond.second_atom_index()]) {
                throw std::invalid_argument{
                    "fixed charge source at molecule " + std::to_string(molecule_index) +
                    ", atom " + std::to_string(bond.second_atom_index()) + " is bonded to atom " +
                    std::to_string(bond.first_atom_index())};
            }
        }
        if (static_cast<std::size_t>(std::ranges::count(selected_mask, true)) ==
            molecule.atom_count()) {
            throw std::invalid_argument{
                "fixed charge embedding leaves no active atoms in molecule " +
                std::to_string(molecule_index)};
        }

        if (molecule.conformer_count() == 0) {
            throw std::invalid_argument{"fixed charge embedding requires a conformer in molecule " +
                                        std::to_string(molecule_index)};
        }
        for (std::size_t conformer_index = 0; conformer_index < molecule.conformer_count();
             ++conformer_index) {
            const auto positions = molecule.conformers()[conformer_index].positions();
            for (std::size_t atom_index = 0; atom_index < positions.size(); ++atom_index) {
                const auto& position = positions[atom_index];
                if (!std::isfinite(position.x) || !std::isfinite(position.y) ||
                    !std::isfinite(position.z)) {
                    throw std::invalid_argument{
                        "fixed charge embedding has non-finite coordinates at molecule " +
                        std::to_string(molecule_index) + ", conformer " +
                        std::to_string(conformer_index) + ", atom " + std::to_string(atom_index)};
                }
            }
            for (auto source = first; source != last; ++source) {
                const auto& source_position = positions[source->atom_index];
                for (std::size_t active_atom = 0; active_atom < positions.size(); ++active_atom) {
                    if (selected_mask[active_atom]) {
                        continue;
                    }
                    const auto& active_position = positions[active_atom];
                    if (source_position.x == active_position.x &&
                        source_position.y == active_position.y &&
                        source_position.z == active_position.z) {
                        throw std::invalid_argument{
                            "fixed charge source at molecule " + std::to_string(molecule_index) +
                            ", conformer " + std::to_string(conformer_index) + ", atom " +
                            std::to_string(source->atom_index) + " coincides with active atom " +
                            std::to_string(active_atom)};
                    }
                }
            }
        }

        auto source_charge_sum = 0.0;
        for (auto source = first; source != last; ++source) {
            source_charge_sum += source->charge;
            if (!std::isfinite(source_charge_sum)) {
                throw std::invalid_argument{"fixed charge source sum is non-finite in molecule " +
                                            std::to_string(molecule_index)};
            }
        }
        target.original_charge = core::total_formal_charge(molecule);
        if (!std::isfinite(target.original_charge)) {
            throw std::invalid_argument{"original formal charge is non-finite in molecule " +
                                        std::to_string(molecule_index)};
        }
        target.active_charge = target.original_charge - source_charge_sum;
        if (!std::isfinite(target.active_charge)) {
            throw std::invalid_argument{"active charge total is non-finite in molecule " +
                                        std::to_string(molecule_index)};
        }
        first = last;
    }

    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        auto& target = targets[molecule_index];
        if (target.sources.empty()) {
            target.original_charge = core::total_formal_charge(molecules[molecule_index]);
            target.active_charge = target.original_charge;
        }
    }
    return targets;
}

} // namespace

auto validate_fixed_charge_embedding(const core::MoleculeCollection& molecules,
                                     const FixedChargeEmbedding& embedding) -> void {
    static_cast<void>(validate_and_group_sources(molecules, embedding));
}

auto make_fixed_charge_partition(const core::MoleculeCollection& molecules,
                                 const FixedChargeEmbedding& embedding) -> FixedChargePartition {
    const auto validated_targets = validate_and_group_sources(molecules, embedding);
    auto active_molecules = std::vector<core::Molecule>{};
    auto targets = std::vector<FixedChargePartitionTarget>(molecules.size());
    active_molecules.reserve(molecules.size());

    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        const auto& source_molecule = molecules[molecule_index];
        const auto& validated = validated_targets[molecule_index];
        auto& target = targets[molecule_index];
        target.sources = validated.sources;
        target.original_charge = validated.original_charge;
        target.active_charge = validated.active_charge;

        if (validated.sources.empty()) {
            target.active_atom_indices.resize(source_molecule.atom_count());
            std::ranges::iota(target.active_atom_indices, std::size_t{0});
            target.active_bond_indices.resize(source_molecule.bond_count());
            std::ranges::iota(target.active_bond_indices, std::size_t{0});
            active_molecules.push_back(source_molecule);
            continue;
        }

        const auto no_active_index = std::numeric_limits<std::size_t>::max();
        auto original_to_active =
            std::vector<std::size_t>(source_molecule.atom_count(), no_active_index);
        auto active_atoms = std::vector<core::Atom>{};
        active_atoms.reserve(source_molecule.atom_count() - validated.sources.size());
        auto source_iterator = validated.sources.begin();
        for (std::size_t original_index = 0; original_index < source_molecule.atom_count();
             ++original_index) {
            if (source_iterator != validated.sources.end() &&
                source_iterator->atom_index == original_index) {
                ++source_iterator;
                continue;
            }
            original_to_active[original_index] = target.active_atom_indices.size();
            target.active_atom_indices.push_back(original_index);
            active_atoms.push_back(source_molecule.atom(original_index));
        }

        auto active_bonds = std::vector<core::Bond>{};
        active_bonds.reserve(source_molecule.bond_count());
        for (std::size_t original_bond_index = 0;
             original_bond_index < source_molecule.bond_count(); ++original_bond_index) {
            const auto& original_bond = source_molecule.bond(original_bond_index);
            const auto active_first = original_to_active[original_bond.first_atom_index()];
            const auto active_second = original_to_active[original_bond.second_atom_index()];
            if (active_first == no_active_index || active_second == no_active_index) {
                throw std::logic_error{
                    "validated fixed-charge partition contains a bond to a fixed source"};
            }
            active_bonds.emplace_back(active_first, active_second, original_bond.order());
            target.active_bond_indices.push_back(original_bond_index);
        }

        auto active_conformers = std::vector<core::Conformer>{};
        active_conformers.reserve(source_molecule.conformer_count());
        target.source_positions.reserve(source_molecule.conformer_count());
        for (const auto& original_conformer : source_molecule.conformers()) {
            auto active_positions = std::vector<core::Position>{};
            active_positions.reserve(target.active_atom_indices.size());
            for (const auto original_atom_index : target.active_atom_indices) {
                active_positions.push_back(original_conformer.positions()[original_atom_index]);
            }
            active_conformers.emplace_back(std::move(active_positions),
                                           std::string{original_conformer.name()});

            auto source_positions = std::vector<core::Position>{};
            source_positions.reserve(validated.sources.size());
            for (const auto& source : validated.sources) {
                source_positions.push_back(original_conformer.positions()[source.atom_index]);
            }
            target.source_positions.push_back(std::move(source_positions));
        }
        active_molecules.emplace_back(std::move(active_atoms), std::move(active_bonds),
                                      std::move(active_conformers),
                                      std::string{source_molecule.name()});
    }

    return {core::MoleculeCollection{std::move(active_molecules), std::string{molecules.name()}},
            std::move(targets), embedding.charge_provenance};
}

} // namespace chargefw::calculation::detail
