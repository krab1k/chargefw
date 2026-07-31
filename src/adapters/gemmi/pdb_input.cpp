#include <chargefw/adapters/gemmi/pdb_input.h>

#include "adapters/native/common_input.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/position.h>

#include <gemmi/pdb.hpp>

#include <cstddef>
#include <istream>
#include <iterator>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::adapters::gemmi::pdb_input {
namespace {

namespace common = chargefw::adapters::native::common_input;

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
                    throw std::runtime_error{"PDB atom '" + atom.name + "' has no known element"};
                }

                result.atoms.emplace_back(atomic_number, selected->charge, atom.name);
                result.positions.push_back(core::Position{
                    .x = selected->pos.x, .y = selected->pos.y, .z = selected->pos.z});
            }
        }
    }

    if (result.atoms.empty()) {
        throw std::runtime_error{"PDB model contains no atoms"};
    }

    return result;
}

auto validate_topology(const std::vector<core::Atom>& reference,
                       const std::vector<core::Atom>& atoms) -> void {
    if (reference.size() != atoms.size()) {
        throw std::runtime_error{"PDB models do not contain the same selected atom sequence"};
    }

    for (std::size_t index = 0; index < reference.size(); ++index) {
        if (reference[index].atomic_number() != atoms[index].atomic_number() ||
            reference[index].formal_charge() != atoms[index].formal_charge() ||
            reference[index].name() != atoms[index].name()) {
            throw std::runtime_error{"PDB models do not contain the same selected atom sequence"};
        }
    }
}

} // namespace

PdbReader::PdbReader(std::istream& input, std::string source, const RecordSelection selection) {
    std::string contents{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    if (input.bad()) {
        throw std::runtime_error{"failed to read PDB input"};
    }

    const auto structure = ::gemmi::read_pdb_string(contents, source);
    if (structure.models.empty()) {
        throw std::runtime_error{"PDB input contains no models"};
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

    const auto name = structure.name.empty() ? source : structure.name;
    record_ = common::make_record(
        std::move(first.atoms), {}, std::move(conformers),
        MoleculeRecordIdentity{.source = std::move(source), .record_index = 0, .record_id = name},
        name);
}

auto PdbReader::next() -> std::optional<ImportedMoleculeRecord> {
    return std::exchange(record_, std::nullopt);
}

} // namespace chargefw::adapters::gemmi::pdb_input
