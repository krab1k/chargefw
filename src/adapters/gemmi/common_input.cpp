#include "common_input.h"

#include "bonds.h"

#include "adapters/native/common_input.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/position.h>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace chargefw::adapters::gemmi::common_input {
namespace {

namespace native_common = chargefw::adapters::native::common_input;

struct ModelAtoms {
    std::vector<core::Atom> atoms;
    std::vector<core::Position> positions;
};

[[nodiscard]] auto include_residue(const ::gemmi::Residue& residue, const RecordSelection selection)
    -> bool {
    switch (selection) {
    case RecordSelection::all:
        return true;
    case RecordSelection::polymers_and_ligands:
        return residue.het_flag != 'H' || residue.name != "HOH";
    case RecordSelection::polymers:
        return residue.het_flag != 'H';
    }

    return false;
}

[[nodiscard]] auto selected_atoms(const ::gemmi::Model& model, const RecordSelection selection)
    -> ModelAtoms {
    ModelAtoms result;

    for (const auto& chain : model.chains) {
        for (const auto& residue : chain.residues) {
            if (!include_residue(residue, selection)) {
                continue;
            }

            for (std::size_t index = 0; index < residue.atoms.size(); ++index) {
                const auto& atom = residue.atoms[index];
                bool previously_seen = false;
                for (std::size_t previous = 0; previous < index; ++previous) {
                    if (residue.atoms[previous].name == atom.name) {
                        previously_seen = true;
                        break;
                    }
                }
                if (previously_seen) {
                    continue;
                }

                const ::gemmi::Atom* selected = std::addressof(atom);
                for (std::size_t candidate_index = index + 1;
                     candidate_index < residue.atoms.size(); ++candidate_index) {
                    const auto& candidate = residue.atoms[candidate_index];
                    if (candidate.name != atom.name) {
                        continue;
                    }
                    if (selected->altloc != '\0' && candidate.altloc == '\0') {
                        selected = std::addressof(candidate);
                    } else if (selected->altloc != '\0' && candidate.altloc == 'A') {
                        selected = std::addressof(candidate);
                    }
                }

                const auto atomic_number = selected->element.atomic_number();
                if (atomic_number <= 0) {
                    throw std::runtime_error{"structural atom '" + atom.name +
                                             "' has no known element"};
                }

                result.atoms.emplace_back(atomic_number, selected->charge, atom.name);
                result.positions.push_back(core::Position{
                    .x = selected->pos.x, .y = selected->pos.y, .z = selected->pos.z});
            }
        }
    }

    if (result.atoms.empty()) {
        throw std::runtime_error{"structural model contains no selected atoms"};
    }

    return result;
}

auto validate_topology(const std::vector<core::Atom>& reference,
                       const std::vector<core::Atom>& atoms) -> void {
    if (reference.size() != atoms.size()) {
        throw std::runtime_error{
            "structural models do not contain the same selected atom sequence"};
    }

    for (std::size_t index = 0; index < reference.size(); ++index) {
        if (reference[index].atomic_number() != atoms[index].atomic_number() ||
            reference[index].formal_charge() != atoms[index].formal_charge() ||
            reference[index].name() != atoms[index].name()) {
            throw std::runtime_error{
                "structural models do not contain the same selected atom sequence"};
        }
    }
}

} // namespace

auto make_record(const ::gemmi::Structure& structure, MoleculeRecordIdentity identity,
                 const RecordSelection selection, const BondStrategy bond_strategy,
                 ::gemmi::cif::Block* const mmcif_block, std::string name)
    -> ImportedMoleculeRecord {
    if (structure.models.empty()) {
        throw std::runtime_error{"structural input contains no models"};
    }

    auto first = selected_atoms(structure.models.front(), selection);
    std::vector<core::Conformer> conformers;
    conformers.reserve(structure.models.size());
    conformers.emplace_back(std::move(first.positions),
                            std::to_string(structure.models.front().num));

    for (std::size_t index = 1; index < structure.models.size(); ++index) {
        auto model = selected_atoms(structure.models[index], selection);
        validate_topology(first.atoms, model.atoms);
        conformers.emplace_back(std::move(model.positions),
                                std::to_string(structure.models[index].num));
    }

    if (name.empty()) {
        name = structure.name;
    }

    std::vector<core::Bond> bonds;
    if (bond_strategy == BondStrategy::templates || bond_strategy == BondStrategy::hybrid) {
        bonds = ::chargefw::adapters::gemmi::bonds::assign_template_bonds(structure.models.front(),
                                                                          selection);
    }
    if (bond_strategy == BondStrategy::explicit_bonds || bond_strategy == BondStrategy::hybrid) {
        const auto explicit_bonds =
            mmcif_block == nullptr
                ? ::chargefw::adapters::gemmi::bonds::assign_explicit_pdb_bonds(structure,
                                                                                selection)
                : ::chargefw::adapters::gemmi::bonds::assign_explicit_mmcif_bonds(
                      structure, *mmcif_block, selection);
        for (const auto& bond : explicit_bonds) {
            bool duplicate = false;
            for (const auto& existing : bonds) {
                if ((existing.first_atom_index() == bond.first_atom_index() &&
                     existing.second_atom_index() == bond.second_atom_index()) ||
                    (existing.first_atom_index() == bond.second_atom_index() &&
                     existing.second_atom_index() == bond.first_atom_index())) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                bonds.push_back(bond);
            }
        }
    }

    return native_common::make_record(std::move(first.atoms), std::move(bonds),
                                      std::move(conformers), std::move(identity), std::move(name));
}

} // namespace chargefw::adapters::gemmi::common_input
