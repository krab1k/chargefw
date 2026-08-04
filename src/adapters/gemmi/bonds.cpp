#include "bonds.h"

#include "component_templates.h"
#include "selection.h"

#include <gemmi/cif.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace chargefw::adapters::gemmi::bonds {
namespace {

class BondAccumulator {
  public:
    void add(const std::size_t first, const std::size_t second, const core::BondOrder order) {
        if (first == second) {
            return;
        }

        const auto key = BondKey{std::min(first, second), std::max(first, second)};
        if (seen_.insert(key).second) {
            bonds_.emplace_back(first, second, order);
        }
    }

    [[nodiscard]] auto take() -> std::vector<core::Bond> {
        return std::move(bonds_);
    }

  private:
    struct BondKey {
        std::size_t first;
        std::size_t second;

        [[nodiscard]] auto operator==(const BondKey&) const -> bool = default;
    };

    struct BondKeyHash {
        [[nodiscard]] auto operator()(const BondKey& key) const noexcept -> std::size_t {
            return std::hash<std::size_t>{}(key.first) ^
                   (std::hash<std::size_t>{}(key.second) << 1U);
        }
    };

    std::vector<core::Bond> bonds_;
    std::unordered_set<BondKey, BondKeyHash> seen_;
};

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

void add_structure_connections(BondAccumulator& bonds, const ::gemmi::Structure& structure,
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
            bonds.add(*first, *second, core::BondOrder::SINGLE);
        }
    }
}

void add_sequential_bonds(BondAccumulator& bonds,
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

        const auto first = previous.find_atom(previous_atom);
        const auto second = current.find_atom(current_atom);
        if (first.has_value() && second.has_value()) {
            bonds.add(*first, *second, core::BondOrder::SINGLE);
        }
    }
}

[[nodiscard]] auto bond_order(const std::string_view value) -> std::optional<core::BondOrder> {
    if (value == "SING") {
        return core::BondOrder::SINGLE;
    }
    if (value == "DOUB") {
        return core::BondOrder::DOUBLE;
    }
    if (value == "TRIP") {
        return core::BondOrder::TRIPLE;
    }
    if (value == "AROM") {
        return core::BondOrder::SINGLE;
    }
    return std::nullopt;
}

void add_component_bonds(BondAccumulator& result, const selection::SelectedResidue& residue,
                         const std::span<const component_templates::BondTemplate> bonds) {
    for (const auto& bond : bonds) {
        const auto first = residue.find_atom(bond.first);
        const auto second = residue.find_atom(bond.second);
        if (first.has_value() && second.has_value()) {
            result.add(*first, *second, bond.order);
        }
    }
}

[[nodiscard]] auto assign_template_bonds(const selection::SelectedModel& model)
    -> std::vector<core::Bond> {
    BondAccumulator result;
    const auto& residues = model.residues();

    for (const auto& residue : residues) {
        const auto component_template = component_templates::find(residue.residue->name);
        if (!component_template) {
            continue;
        }

        add_component_bonds(result, residue, component_template->bonds);
    }

    add_sequential_bonds(result, residues, component_templates::ComponentKind::amino_acid, "C",
                         "N");
    add_sequential_bonds(result, residues, component_templates::ComponentKind::nucleotide, "O3'",
                         "P");

    return result.take();
}

} // namespace

auto explicit_pdb(const ::gemmi::Structure& structure, const selection::SelectedModel& selected)
    -> std::vector<core::Bond> {
    if (structure.models.empty()) {
        return {};
    }

    BondAccumulator result;
    add_structure_connections(result, structure, selected);

    for (const auto& [serial, partners] : structure.conect_map) {
        const auto first = selected.atom_index_by_serial(serial);
        if (!first.has_value()) {
            continue;
        }

        for (const auto partner : partners) {
            const auto second = selected.atom_index_by_serial(partner);
            if (second.has_value()) {
                result.add(*first, *second, core::BondOrder::SINGLE);
            }
        }
    }

    return result.take();
}

auto explicit_mmcif(const ::gemmi::Structure& structure, ::gemmi::cif::Block& block,
                    const selection::SelectedModel& selected) -> std::vector<core::Bond> {
    if (structure.models.empty()) {
        return {};
    }

    const auto& residues = selected.residues();
    BondAccumulator result;

    std::unordered_map<std::string_view, std::vector<component_templates::BondTemplate>>
        component_bonds;
    auto rows =
        block.find("_chem_comp_bond.", {"comp_id", "atom_id_1", "atom_id_2", "value_order"});
    for (const auto row : rows) {
        const std::string_view component{row[0]};
        if (const auto order = bond_order(row[3])) {
            component_bonds[component].emplace_back(component_templates::BondTemplate{
                .first = row[1], .second = row[2], .order = *order});
        }
    }

    for (const auto& residue : residues) {
        const auto found = component_bonds.find(residue.residue->name);
        if (found != component_bonds.end()) {
            add_component_bonds(result, residue, found->second);
        }
    }

    add_structure_connections(result, structure, selected);

    return result.take();
}

auto assign(const selection::SelectedModel& model, const BondStrategy strategy,
            std::vector<core::Bond> explicit_bonds) -> std::vector<core::Bond> {
    switch (strategy) {
    case BondStrategy::none:
        return {};
    case BondStrategy::templates:
        return assign_template_bonds(model);
    case BondStrategy::explicit_bonds:
        return explicit_bonds;
    case BondStrategy::hybrid: {
        BondAccumulator bonds;
        for (const auto& bond : explicit_bonds) {
            bonds.add(bond.first_atom_index(), bond.second_atom_index(), bond.order());
        }

        // Explicit connectivity overrides a conflicting template bond order.
        for (const auto& bond : assign_template_bonds(model)) {
            bonds.add(bond.first_atom_index(), bond.second_atom_index(), bond.order());
        }

        return bonds.take();
    }
    }

    std::abort();
}

} // namespace chargefw::adapters::gemmi::bonds
