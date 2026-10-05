#include "fixed_charge_groups.h"

#include "component.h"

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace chargefw::adapters::detail {
namespace {

struct Ion {
    std::string_view component;
    std::string_view atom;
    int atomic_number;
    double charge;
};

// CCD snapshot records: https://files.rcsb.org/ligands/view/{NA,K,MG,CA,CL,ZN,FE,FE2}.cif
// Verified against components.gz on 2026-10-05. CCD component and atom labels are identical.
constexpr std::array ions{
    Ion{"NA", "NA", 11, 1.0}, Ion{"K", "K", 19, 1.0},     Ion{"MG", "MG", 12, 2.0},
    Ion{"CA", "CA", 20, 2.0}, Ion{"CL", "CL", 17, -1.0},  Ion{"ZN", "ZN", 30, 2.0},
    Ion{"FE", "FE", 26, 3.0}, Ion{"FE2", "FE2", 26, 2.0},
};

constexpr std::string_view provenance = "chargefw:fixed-charge-ions:v1";

auto context(const std::size_t molecule_index, const std::size_t atom_index,
             const std::string_view component, const std::string_view message) -> std::string {
    auto result = std::string{message} + " at record " + std::to_string(molecule_index);
    if (atom_index != std::string_view::npos) {
        result += ", atom " + std::to_string(atom_index);
    }
    if (!component.empty()) {
        result += " (component " + std::string{component} + ")";
    }
    return result;
}

} // namespace

auto resolve_fixed_charge_groups(const std::span<const ImportedMoleculeRecord> records,
                                 const std::span<const std::string> component_ids)
    -> calculation::FixedChargeEmbedding {
    calculation::FixedChargeEmbedding result;
    if (component_ids.empty()) {
        return result;
    }

    std::vector<const Ion*> selected;
    for (const auto& id : component_ids) {
        const auto ion = std::ranges::find(ions, id, &Ion::component);
        if (ion == ions.end()) {
            throw std::invalid_argument("unknown fixed-charge component ID: " + id);
        }
        if (std::ranges::find(selected, std::addressof(*ion)) == selected.end()) {
            selected.push_back(std::addressof(*ion));
        }
    }

    for (std::size_t molecule_index = 0; molecule_index < records.size(); ++molecule_index) {
        const auto& record = records[molecule_index];
        if (record.molecule.atom_count() == 0) {
            continue;
        }
        if (!record.import_metadata ||
            record.import_metadata->atoms.size() != record.molecule.atom_count() ||
            record.import_metadata->components.empty()) {
            throw std::invalid_argument(
                context(molecule_index, std::string_view::npos, {},
                        "fixed-charge resolution requires complete component metadata"));
        }
        const auto& atoms = record.import_metadata->atoms;
        const auto atom_count = record.molecule.atom_count();
        constexpr auto missing_component = std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> component_for_atom(atom_count, missing_component);
        std::vector<bool> bonds(record.molecule.atom_count(), false);
        for (const auto& bond : record.molecule.bonds()) {
            bonds[bond.first_atom_index()] = true;
            bonds[bond.second_atom_index()] = true;
        }
        for (std::size_t component_index = 0;
             component_index < record.import_metadata->components.size(); ++component_index) {
            const auto& component = record.import_metadata->components[component_index];
            if (component.component_id.empty() || component.atom_indices.empty()) {
                throw std::invalid_argument(context(molecule_index, std::string_view::npos,
                                                    component.component_id,
                                                    "component metadata entry is empty"));
            }
            for (const auto atom_index : component.atom_indices) {
                if (atom_index >= atom_count ||
                    component_for_atom[atom_index] != missing_component) {
                    throw std::invalid_argument(
                        context(molecule_index, atom_index, component.component_id,
                                "component metadata has an invalid or repeated atom index"));
                }
                component_for_atom[atom_index] = component_index;
            }
        }
        if (std::ranges::find(component_for_atom, missing_component) != component_for_atom.end()) {
            throw std::invalid_argument(context(molecule_index, std::string_view::npos, {},
                                                "component metadata does not cover every atom"));
        }

        const auto first_source = result.sources.size();
        for (const auto& component : record.import_metadata->components) {
            const auto canonical_name = canonical_component_name(
                atoms[component.atom_indices.front()].structural_labels, component.component_id);
            const auto target = std::ranges::find(selected, canonical_name,
                                                  [](const Ion* ion) { return ion->component; });
            if (target == selected.end()) {
                continue;
            }
            if (component.atom_indices.size() != 1) {
                throw std::invalid_argument(
                    context(molecule_index, component.atom_indices.front(), component.component_id,
                            "selected fixed-charge component is not monatomic"));
            }
            auto component_atoms = std::vector<ComponentAtom>{};
            component_atoms.reserve(component.atom_indices.size());
            for (const auto atom_index : component.atom_indices) {
                const auto source_atom_name = canonical_atom_name(
                    atoms[atom_index].structural_labels, record.molecule.atom(atom_index).name());
                component_atoms.push_back({std::string{source_atom_name}, atom_index});
            }
            const auto view = ComponentView{canonical_name, component_atoms};
            const auto atom_index = view.find_atom((*target)->atom);
            if (!atom_index.has_value()) {
                throw std::invalid_argument(
                    context(molecule_index, component.atom_indices.front(), component.component_id,
                            "selected fixed-charge component atom identity mismatch"));
            }
            const auto& atom = record.molecule.atom(*atom_index);
            if (atom.atomic_number() != (*target)->atomic_number) {
                throw std::invalid_argument(
                    context(molecule_index, *atom_index, component.component_id,
                            "selected fixed-charge component atom identity mismatch"));
            }
            if (bonds[*atom_index]) {
                throw std::invalid_argument(
                    context(molecule_index, *atom_index, component.component_id,
                            "selected fixed-charge component has an incident graph bond"));
            }
            result.sources.push_back({molecule_index, *atom_index, (*target)->charge});
        }
        std::ranges::sort(result.sources.begin() + static_cast<std::ptrdiff_t>(first_source),
                          result.sources.end(), {}, &calculation::FixedAtomCharge::atom_index);
    }
    if (!result.sources.empty()) {
        result.charge_provenance = provenance;
    }
    return result;
}

} // namespace chargefw::adapters::detail
