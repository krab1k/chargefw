#include "structure_import.h"

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

namespace chargefw::adapters::gemmi::structure_import {
namespace {

namespace native_common = chargefw::adapters::native::common_input;

struct ModelAtoms {
    std::vector<core::Atom> atoms;
    std::vector<core::Position> positions;
};

[[nodiscard]] auto model_atoms(const selection::SelectedModel& model) -> ModelAtoms {
    ModelAtoms result;
    result.atoms.reserve(model.atoms().size());
    result.positions.reserve(model.atoms().size());

    for (const auto& selected_atom : model.atoms()) {
        const auto& atom = *selected_atom.atom;
        const auto atomic_number = atom.element.atomic_number();
        if (atomic_number <= 0) {
            throw std::runtime_error{"structural atom '" + atom.name + "' has no known element"};
        }

        result.atoms.emplace_back(atomic_number, atom.charge, atom.name);
        result.positions.push_back(
            core::Position{.x = atom.pos.x, .y = atom.pos.y, .z = atom.pos.z});
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
    auto first = model_atoms(selected_model);
    std::vector<core::Conformer> conformers;
    conformers.reserve(structure.models.size());
    conformers.emplace_back(std::move(first.positions),
                            std::to_string(structure.models.front().num));

    for (std::size_t index = 1; index < structure.models.size(); ++index) {
        const auto selected = ::chargefw::adapters::gemmi::selection::SelectedModel{
            structure.models[index], selection};
        auto model = model_atoms(selected);
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

} // namespace chargefw::adapters::gemmi::structure_import
