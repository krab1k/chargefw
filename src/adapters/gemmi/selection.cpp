#include "selection.h"

namespace chargefw::adapters::gemmi::selection {

auto include_residue(const ::gemmi::Residue& residue, const RecordSelection selection) -> bool {
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

auto is_first_named_atom(const ::gemmi::Residue& residue, const std::size_t atom_index) -> bool {
    const auto& atom = residue.atoms[atom_index];
    for (std::size_t previous = 0; previous < atom_index; ++previous) {
        if (residue.atoms[previous].name == atom.name) {
            return false;
        }
    }
    return true;
}

auto select_altloc(const ::gemmi::Residue& residue, const std::size_t first_atom_index)
    -> const ::gemmi::Atom& {
    const auto& first = residue.atoms[first_atom_index];
    const auto* selected = std::addressof(first);

    for (std::size_t candidate_index = first_atom_index + 1; candidate_index < residue.atoms.size();
         ++candidate_index) {
        const auto& candidate = residue.atoms[candidate_index];
        if (candidate.name != first.name) {
            continue;
        }

        if (selected->altloc != '\0' && candidate.altloc == '\0') {
            selected = std::addressof(candidate);
        } else if (selected->altloc != '\0' && selected->altloc != 'A' && candidate.altloc == 'A') {
            selected = std::addressof(candidate);
        }
    }

    return *selected;
}

SelectedModel::SelectedModel(const ::gemmi::Model& model, const RecordSelection selection) {
    std::size_t atom_index = 0;
    for (const auto& chain : model.chains) {
        for (const auto& residue : chain.residues) {
            if (!include_residue(residue, selection)) {
                continue;
            }

            SelectedResidue selected{.residue = std::addressof(residue), .atom_indices = {}};
            for (std::size_t index = 0; index < residue.atoms.size(); ++index) {
                if (!is_first_named_atom(residue, index)) {
                    continue;
                }

                const auto& atom = select_altloc(residue, index);
                atoms_.push_back({.atom = std::addressof(atom), .index = atom_index});
                selected.atom_indices.emplace_back(atom.name, atom_index++);
            }
            residues_.push_back(std::move(selected));
        }
    }
}

auto SelectedModel::atoms() const -> const std::vector<SelectedAtom>& {
    return atoms_;
}

auto SelectedModel::residues() const -> const std::vector<SelectedResidue>& {
    return residues_;
}

auto SelectedModel::atom_index(const ::gemmi::Atom* atom) const -> std::optional<std::size_t> {
    for (const auto& selected : atoms_) {
        if (selected.atom == atom) {
            return selected.index;
        }
    }
    return std::nullopt;
}

auto atom_index(const SelectedResidue& residue, const std::string_view name)
    -> std::optional<std::size_t> {
    for (const auto& [atom_name, index] : residue.atom_indices) {
        if (atom_name == name) {
            return index;
        }
    }
    return std::nullopt;
}

} // namespace chargefw::adapters::gemmi::selection
