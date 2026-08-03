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

} // namespace chargefw::adapters::gemmi::selection
