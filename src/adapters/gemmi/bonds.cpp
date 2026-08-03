#include "bonds.h"

#include "component_templates.h"
#include "selection.h"

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

[[nodiscard]] auto selected_atom_index(const ::gemmi::Model& model,
                                       const selection::SelectedModel& selected,
                                       const ::gemmi::AtomAddress& address)
    -> std::optional<std::size_t> {
    const auto cra = model.find_cra(address);
    if (cra.atom != nullptr) {
        return selected.atom_index(cra.atom);
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
                    return selected.atom_index(std::addressof(atom));
                }
            }
        }
    }
    return std::nullopt;
}

void add_structure_connections(std::vector<core::Bond>& bonds, const ::gemmi::Structure& structure,
                               const selection::SelectedModel& selected) {
    if (structure.models.empty()) {
        return;
    }

    const auto& model = structure.models.front();
    for (const auto& connection : structure.connections) {
        if (connection.type != ::gemmi::Connection::Covale &&
            connection.type != ::gemmi::Connection::Disulf) {
            continue;
        }

        const auto first = selected_atom_index(model, selected, connection.partner1);
        const auto second = selected_atom_index(model, selected, connection.partner2);
        if (first.has_value() && second.has_value()) {
            add_bond(bonds, *first, *second, core::BondOrder::SINGLE);
        }
    }
}

void add_sequential_bonds(std::vector<core::Bond>& bonds,
                          const std::span<const selection::SelectedResidue> residues,
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

        const auto first = selection::atom_index(previous, previous_atom);
        const auto second = selection::atom_index(current, current_atom);
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

[[nodiscard]] auto assign_template_bonds(const selection::SelectedModel& model)
    -> std::vector<core::Bond> {
    std::vector<core::Bond> result;
    const auto& residues = model.residues();

    for (const auto& residue : residues) {
        const auto component_template = component_templates::find(residue.residue->name);
        if (!component_template) {
            continue;
        }

        for (const auto& template_bond : component_template->bonds) {
            const auto first = selection::atom_index(residue, template_bond.first);
            const auto second = selection::atom_index(residue, template_bond.second);
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

} // namespace

auto explicit_pdb(const ::gemmi::Structure& structure, const selection::SelectedModel& selected)
    -> std::vector<core::Bond> {
    if (structure.models.empty()) {
        return {};
    }

    std::vector<core::Bond> result;
    add_structure_connections(result, structure, selected);

    for (const auto& [serial, partners] : structure.conect_map) {
        const auto first = selected.atom_index_by_serial(serial);
        if (!first.has_value()) {
            continue;
        }

        for (const auto partner : partners) {
            const auto second = selected.atom_index_by_serial(partner);
            if (second.has_value()) {
                add_bond(result, *first, *second, core::BondOrder::SINGLE);
            }
        }
    }

    return result;
}

auto explicit_mmcif(const ::gemmi::Structure& structure, ::gemmi::cif::Block& block,
                    const selection::SelectedModel& selected) -> std::vector<core::Bond> {
    if (structure.models.empty()) {
        return {};
    }

    const auto& residues = selected.residues();
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
            const auto first = selection::atom_index(residue, first_name);
            const auto second = selection::atom_index(residue, second_name);
            if (first.has_value() && second.has_value()) {
                add_bond(result, *first, *second, order);
            }
        }
    }

    add_structure_connections(result, structure, selected);

    return result;
}

auto assign(const selection::SelectedModel& model, const BondStrategy strategy,
            std::vector<core::Bond> explicit_bonds) -> std::vector<core::Bond> {
    if (strategy == BondStrategy::none) {
        return {};
    }

    std::vector<core::Bond> result;
    if (strategy == BondStrategy::templates || strategy == BondStrategy::hybrid) {
        result = assign_template_bonds(model);
    }

    if (strategy == BondStrategy::explicit_bonds || strategy == BondStrategy::hybrid) {
        for (const auto& bond : explicit_bonds) {
            add_bond(result, bond.first_atom_index(), bond.second_atom_index(), bond.order());
        }
    }

    return result;
}

} // namespace chargefw::adapters::gemmi::bonds
