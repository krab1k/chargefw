#include "bonds.h"

#include "component_templates.h"

#include <gemmi/cif.hpp>

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

struct SelectedAtom {
    const ::gemmi::Atom* atom;
    std::size_t index;
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

[[nodiscard]] auto selected_atoms(const ::gemmi::Model& model, const RecordSelection selection)
    -> std::vector<SelectedAtom> {
    std::vector<SelectedAtom> result;
    std::size_t atom_index = 0;

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

                result.push_back({.atom = chosen, .index = atom_index++});
            }
        }
    }

    return result;
}

[[nodiscard]] auto selected_atom_index(const std::vector<SelectedAtom>& atoms,
                                       const ::gemmi::Atom* atom) -> std::optional<std::size_t> {
    for (const auto& selected : atoms) {
        if (selected.atom == atom) {
            return selected.index;
        }
    }
    return std::nullopt;
}

[[nodiscard]] auto selected_atom_index(const ::gemmi::Model& model,
                                       const std::vector<SelectedAtom>& atoms,
                                       const ::gemmi::AtomAddress& address)
    -> std::optional<std::size_t> {
    const auto cra = model.find_cra(address);
    if (cra.atom != nullptr) {
        return selected_atom_index(atoms, cra.atom);
    }

    for (const auto& chain : model.chains) {
        if (chain.name != address.chain_name) {
            continue;
        }
        for (const auto& residue : chain.residues) {
            if (residue.seqid != address.res_id.seqid ||
                (!address.res_id.name.empty() && residue.name != address.res_id.name)) {
                continue;
            }
            for (const auto& atom : residue.atoms) {
                if (atom.name == address.atom_name &&
                    (address.altloc == '\0' || atom.altloc == address.altloc)) {
                    return selected_atom_index(atoms, std::addressof(atom));
                }
            }
        }
    }
    return std::nullopt;
}

[[nodiscard]] auto atom_by_serial(const ::gemmi::Model& model, const int serial)
    -> const ::gemmi::Atom* {
    for (const auto& chain : model.chains) {
        for (const auto& residue : chain.residues) {
            for (const auto& atom : residue.atoms) {
                if (atom.serial == serial) {
                    return std::addressof(atom);
                }
            }
        }
    }
    return nullptr;
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

[[nodiscard]] auto bond_order(const std::string_view value) -> core::BondOrder {
    if (value == "SING") {
        return core::BondOrder::SINGLE;
    }
    if (value == "DOUB") {
        return core::BondOrder::DOUBLE;
    }
    if (value == "TRIP") {
        return core::BondOrder::TRIPLE;
    }
    return core::BondOrder::UNKNOWN;
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

auto assign_explicit_pdb_bonds(const ::gemmi::Structure& structure, const RecordSelection selection)
    -> std::vector<core::Bond> {
    if (structure.models.empty()) {
        return {};
    }

    const auto& model = structure.models.front();
    const auto atoms = selected_atoms(model, selection);
    std::vector<core::Bond> result;

    for (const auto& connection : structure.connections) {
        if (connection.type != ::gemmi::Connection::Covale &&
            connection.type != ::gemmi::Connection::Disulf) {
            continue;
        }

        const auto first = selected_atom_index(model, atoms, connection.partner1);
        const auto second = selected_atom_index(model, atoms, connection.partner2);
        if (first.has_value() && second.has_value()) {
            add_bond(result, *first, *second, core::BondOrder::SINGLE);
        }
    }

    for (const auto& [serial, partners] : structure.conect_map) {
        const auto* first_atom = atom_by_serial(model, serial);
        if (first_atom == nullptr) {
            continue;
        }

        const auto first = selected_atom_index(atoms, first_atom);
        if (!first.has_value()) {
            continue;
        }

        for (const auto partner : partners) {
            const auto* second_atom = atom_by_serial(model, partner);
            const auto second =
                second_atom == nullptr ? std::nullopt : selected_atom_index(atoms, second_atom);
            if (second.has_value()) {
                add_bond(result, *first, *second, core::BondOrder::SINGLE);
            }
        }
    }

    return result;
}

auto assign_explicit_mmcif_bonds(const ::gemmi::Structure& structure, ::gemmi::cif::Block& block,
                                 const RecordSelection selection) -> std::vector<core::Bond> {
    if (structure.models.empty()) {
        return {};
    }

    const auto& model = structure.models.front();
    const auto residues = selected_residue_atoms(model, selection);
    std::vector<core::Bond> result;

    auto component_bonds =
        block.find("_chem_comp_bond.", {"comp_id", "atom_id_1", "atom_id_2", "value_order"});
    for (const auto row : component_bonds) {
        const std::string_view component{row[0]};
        const std::string_view first_name{row[1]};
        const std::string_view second_name{row[2]};
        const auto order = bond_order(row[3]);
        for (const auto& residue : residues) {
            if (residue.residue->name != component) {
                continue;
            }
            const auto first = atom_index(residue, first_name);
            const auto second = atom_index(residue, second_name);
            if (first.has_value() && second.has_value()) {
                add_bond(result, *first, *second, order);
            }
        }
    }

    const auto atoms = selected_atoms(model, selection);
    for (const auto& connection : structure.connections) {
        if (connection.type != ::gemmi::Connection::Covale &&
            connection.type != ::gemmi::Connection::Disulf) {
            continue;
        }

        const auto first = selected_atom_index(model, atoms, connection.partner1);
        const auto second = selected_atom_index(model, atoms, connection.partner2);
        if (first.has_value() && second.has_value()) {
            add_bond(result, *first, *second, core::BondOrder::SINGLE);
        }
    }

    return result;
}

} // namespace chargefw::adapters::gemmi::bonds
