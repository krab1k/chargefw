#include "bonds.h"

#include "component_templates.h"

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::adapters::gemmi::bonds {
namespace {

struct ResidueAtoms {
    const ::gemmi::Residue* residue;
    std::vector<std::pair<std::string, std::size_t>> atom_indices;
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

void add_bond(std::vector<core::Bond>& bonds, const std::size_t first, const std::size_t second,
              const core::BondOrder order) {
    if (first == second) {
        return;
    }

    for (const auto& bond : bonds) {
        if ((bond.first_atom_index() == first && bond.second_atom_index() == second) ||
            (bond.first_atom_index() == second && bond.second_atom_index() == first)) {
            return;
        }
    }

    bonds.emplace_back(first, second, order);
}

[[nodiscard]] auto selected_residue_atoms(const ::gemmi::Model& model,
                                          const RecordSelection selection)
    -> std::vector<ResidueAtoms> {
    std::vector<ResidueAtoms> result;
    std::size_t atom_index = 0;

    for (const auto& chain : model.chains) {
        for (const auto& residue : chain.residues) {
            if (!include_residue(residue, selection)) {
                continue;
            }

            ResidueAtoms selected{.residue = std::addressof(residue), .atom_indices = {}};
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

                const ::gemmi::Atom* chosen = std::addressof(atom);
                for (std::size_t candidate_index = index + 1;
                     candidate_index < residue.atoms.size(); ++candidate_index) {
                    const auto& candidate = residue.atoms[candidate_index];
                    if (candidate.name != atom.name) {
                        continue;
                    }
                    if (chosen->altloc != '\0' && candidate.altloc == '\0') {
                        chosen = std::addressof(candidate);
                    } else if (chosen->altloc != '\0' && candidate.altloc == 'A') {
                        chosen = std::addressof(candidate);
                    }
                }

                selected.atom_indices.emplace_back(chosen->name, atom_index++);
            }
            result.push_back(std::move(selected));
        }
    }

    return result;
}

[[nodiscard]] auto atom_index(const ResidueAtoms& residue, const std::string_view name)
    -> std::optional<std::size_t> {
    for (const auto& [atom_name, index] : residue.atom_indices) {
        if (atom_name == name) {
            return index;
        }
    }
    return std::nullopt;
}

void add_sequential_bonds(std::vector<core::Bond>& bonds,
                          const std::span<const ResidueAtoms> residues,
                          const component_templates::ComponentKind kind,
                          const std::string_view previous_atom,
                          const std::string_view current_atom) {
    for (std::size_t index = 1; index < residues.size(); ++index) {
        const auto& previous = residues[index - 1];
        const auto& current = residues[index];
        const auto previous_template = component_templates::find(previous.residue->name);
        const auto current_template = component_templates::find(current.residue->name);
        if (!previous_template || !current_template || previous_template->kind != kind ||
            current_template->kind != kind ||
            previous.residue->seqid.num.value + 1 != current.residue->seqid.num.value) {
            continue;
        }

        const auto first = atom_index(previous, previous_atom);
        const auto second = atom_index(current, current_atom);
        if (first.has_value() && second.has_value()) {
            add_bond(bonds, *first, *second, core::BondOrder::SINGLE);
        }
    }
}

} // namespace

auto assign_template_bonds(const ::gemmi::Model& model, const RecordSelection selection)
    -> std::vector<core::Bond> {
    std::vector<core::Bond> result;
    const auto residues = selected_residue_atoms(model, selection);

    for (const auto& residue : residues) {
        const auto component_template = component_templates::find(residue.residue->name);
        if (!component_template) {
            continue;
        }

        for (const auto& template_bond : component_template->bonds) {
            const auto first = atom_index(residue, template_bond.first);
            const auto second = atom_index(residue, template_bond.second);
            if (first.has_value() && second.has_value()) {
                add_bond(result, *first, *second, template_bond.order);
            }
        }
    }

    add_sequential_bonds(result, residues, component_templates::ComponentKind::amino_acid, "C",
                         "N");
    add_sequential_bonds(result, residues, component_templates::ComponentKind::nucleotide, "O3'",
                         "P");

    return result;
}

} // namespace chargefw::adapters::gemmi::bonds
