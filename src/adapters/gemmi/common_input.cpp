#include "common_input.h"

#include "bonds.h"
#include "selection.h"

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

[[nodiscard]] auto selected_atoms(const ::gemmi::Model& model, const RecordSelection selection)
    -> ModelAtoms {
    ModelAtoms result;

    for (const auto& chain : model.chains) {
        for (const auto& residue : chain.residues) {
            if (!::chargefw::adapters::gemmi::selection::include_residue(residue, selection)) {
                continue;
            }

            for (std::size_t index = 0; index < residue.atoms.size(); ++index) {
                if (!::chargefw::adapters::gemmi::selection::is_first_named_atom(residue, index)) {
                    continue;
                }

                const auto& selected =
                    ::chargefw::adapters::gemmi::selection::select_altloc(residue, index);

                const auto atomic_number = selected.element.atomic_number();
                if (atomic_number <= 0) {
                    throw std::runtime_error{"structural atom '" + selected.name +
                                             "' has no known element"};
                }

                result.atoms.emplace_back(atomic_number, selected.charge, selected.name);
                result.positions.push_back(
                    core::Position{.x = selected.pos.x, .y = selected.pos.y, .z = selected.pos.z});
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

    const auto selected_model =
        ::chargefw::adapters::gemmi::selection::SelectedModel{structure.models.front(), selection};
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

    auto bonds = ::chargefw::adapters::gemmi::bonds::assign(structure, selected_model,
                                                            bond_strategy, mmcif_block);

    return native_common::make_record(std::move(first.atoms), std::move(bonds),
                                      std::move(conformers), std::move(identity), std::move(name));
}

} // namespace chargefw::adapters::gemmi::common_input
