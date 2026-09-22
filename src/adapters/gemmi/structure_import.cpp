#include "structure_import.h"

#include "bonds.h"
#include "selection.h"

#include "adapters/native/common_input.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/position.h>

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <span>
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

[[nodiscard]] auto source_order(const std::span<const SourceAtomReference> references)
    -> std::vector<std::size_t> {
    auto result = std::vector<std::size_t>(references.size());
    std::ranges::iota(result, std::size_t{0});
    std::ranges::sort(result, {}, [&](const auto index) { return references[index].position; });
    for (std::size_t index = 1; index < result.size(); ++index) {
        if (references[result[index - 1]].position == references[result[index]].position) {
            throw std::runtime_error{"structural source mapping contains duplicate atom positions"};
        }
    }
    return result;
}

template <typename T>
[[nodiscard]] auto reordered(std::vector<T> values, const std::span<const std::size_t> order)
    -> std::vector<T> {
    if (values.size() != order.size()) {
        throw std::runtime_error{"structural source mapping size does not match selected atoms"};
    }
    auto result = std::vector<T>{};
    result.reserve(values.size());
    for (const auto index : order) {
        if (index >= values.size()) {
            throw std::runtime_error{"structural source mapping contains an invalid atom index"};
        }
        result.push_back(std::move(values[index]));
    }
    return result;
}

[[nodiscard]] auto import_reference(const selection::SelectedModel& model) -> ModelAtoms {
    ModelAtoms result;
    result.atoms.reserve(model.atoms().size());
    result.positions.reserve(model.atoms().size());

    for (const auto* atom : model.atoms()) {
        const auto atomic_number = atom->element.atomic_number();
        if (atomic_number <= 0) {
            throw std::runtime_error{"structural atom '" + atom->name + "' has no known element"};
        }

        result.atoms.emplace_back(atomic_number, atom->charge, atom->name);
        result.positions.push_back(
            core::Position{.x = atom->pos.x, .y = atom->pos.y, .z = atom->pos.z});
    }

    if (result.atoms.empty()) {
        throw std::runtime_error{"structural model contains no selected atoms"};
    }

    return result;
}

[[nodiscard]] auto conformer_positions(const selection::SelectedModel& model,
                                       const std::span<const core::Atom> reference,
                                       const std::span<const SourceAtomReference> reference_mapping,
                                       const std::span<const SourceAtomReference> source_mapping,
                                       const std::span<const std::size_t> order)
    -> std::vector<core::Position> {
    if (reference.size() != model.atoms().size() || reference.size() != reference_mapping.size() ||
        reference.size() != source_mapping.size() || reference.size() != order.size()) {
        throw std::runtime_error{
            "structural models do not contain the same selected atom sequence"};
    }

    std::vector<core::Position> positions;
    positions.reserve(reference.size());
    for (std::size_t index = 0; index < reference.size(); ++index) {
        const auto source_index = order[index];
        if (source_index >= model.atoms().size()) {
            throw std::runtime_error{"structural source mapping contains an invalid atom index"};
        }
        const auto& source = *model.atoms()[source_index];
        if (reference[index].atomic_number() != source.element.atomic_number() ||
            reference[index].formal_charge() != source.charge ||
            reference[index].name() != source.name) {
            throw std::runtime_error{
                "structural models do not contain the same selected atom sequence"};
        }

        const auto& reference_labels = reference_mapping[index].structural_labels;
        const auto& source_labels = source_mapping[index].structural_labels;
        const auto same_hierarchy = [&] {
            if (!reference_labels.has_value() || !source_labels.has_value()) {
                return reference_labels.has_value() == source_labels.has_value();
            }
            const auto& first = *reference_labels;
            const auto& second = *source_labels;
            return first.author.atom == second.author.atom &&
                   first.author.residue == second.author.residue &&
                   first.author.chain == second.author.chain &&
                   first.author.sequence == second.author.sequence &&
                   first.label.atom == second.label.atom &&
                   first.label.residue == second.label.residue &&
                   first.label.chain == second.label.chain &&
                   first.label.sequence == second.label.sequence && first.entity == second.entity &&
                   first.insertion_code == second.insertion_code && first.segment == second.segment;
        }();
        if (!same_hierarchy) {
            throw std::runtime_error{
                "structural models do not contain the same selected atom identity at atom " +
                std::to_string(index)};
        }

        positions.push_back(
            core::Position{.x = source.pos.x, .y = source.pos.y, .z = source.pos.z});
    }

    return positions;
}

[[nodiscard]] auto remap_bonds(std::vector<core::Bond> bonds,
                               const std::span<const std::size_t> order)
    -> std::vector<core::Bond> {
    auto old_to_new = std::vector<std::size_t>(order.size());
    for (std::size_t new_index = 0; new_index < order.size(); ++new_index) {
        old_to_new[order[new_index]] = new_index;
    }
    for (auto& bond : bonds) {
        if (bond.first_atom_index() >= old_to_new.size() ||
            bond.second_atom_index() >= old_to_new.size()) {
            throw std::runtime_error{"structural bond endpoint is outside selected atom mapping"};
        }
        bond = core::Bond{old_to_new[bond.first_atom_index()], old_to_new[bond.second_atom_index()],
                          bond.order()};
    }
    return bonds;
}

} // namespace

auto make_record(const ::gemmi::Structure& structure,
                 const std::span<const selection::SelectedModel> selected_models,
                 MoleculeRecordIdentity identity, const RecordSelection selection,
                 const BondStrategy bond_strategy, const ConformerSelection conformer_selection,
                 std::vector<core::Bond> explicit_bonds, std::string name,
                 const MolecularSourceFormat format, std::vector<SourceModelMapping> source_models,
                 const SourceConnectivity source_connectivity) -> ImportedMoleculeRecord {
    if (structure.models.empty()) {
        throw std::runtime_error{"structural input contains no models"};
    }

    const auto retained_model_count =
        conformer_selection == ConformerSelection::all ? structure.models.size() : std::size_t{1};
    if (selected_models.size() != retained_model_count) {
        throw std::runtime_error{"structural selection does not match retained models"};
    }

    const auto& selected_model = selected_models.front();
    auto first = import_reference(selected_model);
    if (source_models.size() != retained_model_count || source_models.empty() ||
        source_models.front().conformer.sites.size() != first.atoms.size()) {
        throw std::runtime_error{"structural source mapping does not match selected models"};
    }

    auto source_orders = std::vector<std::vector<std::size_t>>{};
    source_orders.reserve(source_models.size());
    for (auto& source_model : source_models) {
        auto order = source_order(source_model.conformer.sites);
        source_model.conformer.sites = reordered(std::move(source_model.conformer.sites), order);
        source_orders.push_back(std::move(order));
    }
    first.atoms = reordered(std::move(first.atoms), source_orders.front());
    first.positions = reordered(std::move(first.positions), source_orders.front());

    std::vector<core::Conformer> conformers;
    conformers.reserve(retained_model_count);
    conformers.emplace_back(std::move(first.positions),
                            std::to_string(structure.models.front().num));

    for (std::size_t index = 1;
         index < structure.models.size() && conformer_selection == ConformerSelection::all;
         ++index) {
        auto positions = conformer_positions(
            selected_models[index], first.atoms, source_models.front().conformer.sites,
            source_models[index].conformer.sites, source_orders[index]);
        conformers.emplace_back(std::move(positions), std::to_string(structure.models[index].num));
    }

    if (name.empty()) {
        name = structure.name;
    }

    auto bonds = remap_bonds(::chargefw::adapters::gemmi::bonds::assign(
                                 selected_model, bond_strategy, std::move(explicit_bonds)),
                             source_orders.front());

    auto atom_references = source_models.front().conformer.sites;
    auto conformer_references = std::vector<SourceConformerReference>{};
    conformer_references.reserve(source_models.size());
    for (auto& source_model : source_models) {
        conformer_references.push_back(std::move(source_model.conformer));
    }
    auto metadata = MoleculeImportMetadata{
        .format = format,
        .atoms = std::move(atom_references),
        .conformers = std::move(conformer_references),
        .record_selection = std::string{::chargefw::adapters::gemmi::to_string(selection)},
        .alternate_location_selection = "blank-then-A-then-first",
        .conformer_selection = std::string{::chargefw::adapters::to_string(conformer_selection)},
        .bond_strategy = std::string{::chargefw::adapters::gemmi::to_string(bond_strategy)},
        .source_connectivity = source_connectivity,
    };

    return native_common::make_record(std::move(first.atoms), std::move(bonds),
                                      std::move(conformers), std::move(identity), std::move(name),
                                      {}, std::move(metadata));
}

} // namespace chargefw::adapters::gemmi::structure_import
