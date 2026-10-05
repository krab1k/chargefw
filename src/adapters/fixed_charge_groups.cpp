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
    bool common = false;
};

// CCD snapshot records: https://files.rcsb.org/pub/pdb/data/monomers/components.cif.gz
// Verified on 2026-10-05: released monatomic ions with nonzero atom/component charge agreement.
constexpr std::array ions{
    Ion{"NA", "NA", 11, 1.0, true}, Ion{"K", "K", 19, 1.0, true},    Ion{"MG", "MG", 12, 2.0, true},
    Ion{"CA", "CA", 20, 2.0, true}, Ion{"CL", "CL", 17, -1.0, true}, Ion{"ZN", "ZN", 30, 2.0, true},
    Ion{"FE", "FE", 26, 3.0, true}, Ion{"FE2", "FE", 26, 2.0, true}, Ion{"0BE", "BE", 4, 2.0},
    Ion{"3CO", "CO", 27, 3.0},      Ion{"3NI", "NI", 28, 3.0},       Ion{"4MO", "MO", 42, 4.0},
    Ion{"4PU", "PU", 94, 4.0},      Ion{"4TI", "TI", 22, 4.0},       Ion{"6MO", "MO", 42, 6.0},
    Ion{"AG", "AG", 47, 1.0},       Ion{"AL", "AL", 13, 3.0},        Ion{"AM", "AM", 95, 3.0},
    Ion{"AU", "AU", 79, 1.0},       Ion{"AU3", "AU", 79, 3.0},       Ion{"BA", "BA", 56, 2.0},
    Ion{"BR", "BR", 35, -1.0},      Ion{"BS3", "BI", 83, 3.0},       Ion{"CD", "CD", 48, 2.0},
    Ion{"CE", "CE", 58, 3.0},       Ion{"CF", "CF", 98, 3.0},        Ion{"CO", "CO", 27, 2.0},
    Ion{"CR", "CR", 24, 3.0},       Ion{"CS", "CS", 55, 1.0},        Ion{"CU", "CU", 29, 2.0},
    Ion{"CU1", "CU", 29, 1.0},      Ion{"CU3", "CU", 29, 3.0},       Ion{"D8U", "D", 1, 1.0},
    Ion{"DY", "DY", 66, 3.0},       Ion{"ER3", "ER", 68, 3.0},       Ion{"EU", "EU", 63, 2.0},
    Ion{"EU3", "EU", 63, 3.0},      Ion{"F", "F", 9, -1.0},          Ion{"GA", "GA", 31, 3.0},
    Ion{"GD3", "GD", 64, 3.0},      Ion{"HG", "HG", 80, 2.0},        Ion{"HO3", "HO", 67, 3.0},
    Ion{"IN", "IN", 49, 3.0},       Ion{"IOD", "I", 53, -1.0},       Ion{"IR", "IR", 77, 4.0},
    Ion{"IR3", "IR", 77, 3.0},      Ion{"LA", "LA", 57, 3.0},        Ion{"LI", "LI", 3, 1.0},
    Ion{"LU", "LU", 71, 3.0},       Ion{"MN", "MN", 25, 2.0},        Ion{"MN3", "MN", 25, 3.0},
    Ion{"ND", "ND", 60, 3.0},       Ion{"NI", "NI", 28, 2.0},        Ion{"OS", "OS", 76, 3.0},
    Ion{"OS4", "OS", 76, 4.0},      Ion{"PB", "PB", 82, 2.0},        Ion{"PD", "PD", 46, 2.0},
    Ion{"PR", "PR", 59, 3.0},       Ion{"PT", "PT", 78, 2.0},        Ion{"PT4", "PT", 78, 4.0},
    Ion{"RB", "RB", 37, 1.0},       Ion{"RH", "RH1", 45, 1.0},       Ion{"RH3", "RH", 45, 3.0},
    Ion{"RHF", "RH", 45, 2.0},      Ion{"RU", "RU", 44, 3.0},        Ion{"SB", "SB", 51, 3.0},
    Ion{"SM", "SM", 62, 3.0},       Ion{"SR", "SR", 38, 2.0},        Ion{"TB", "TB", 65, 3.0},
    Ion{"TH", "TH", 90, 4.0},       Ion{"TL", "TL", 81, 1.0},        Ion{"V", "V", 23, 3.0},
    Ion{"W", "W", 74, 6.0},         Ion{"Y1", "Y", 39, 2.0},         Ion{"YB", "YB", 70, 3.0},
    Ion{"YB2", "YB", 70, 2.0},      Ion{"YT3", "Y", 39, 3.0},        Ion{"ZCM", "CM", 96, 3.0},
    Ion{"ZR", "ZR", 40, 4.0},       Ion{"ZTM", "AC", 89, 3.0},
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

auto fixed_charge_ion_names(const bool common_only) -> std::vector<std::string> {
    std::vector<std::string> names;
    for (const auto& ion : ions) {
        if (!common_only || ion.common) {
            names.emplace_back(ion.component);
        }
    }
    if (!common_only) {
        std::ranges::sort(names);
    }
    return names;
}

auto resolve_fixed_charge_groups(const std::span<const ImportedMoleculeRecord> records,
                                 const std::span<const std::string> component_ids)
    -> calculation::FixedChargeGroups {
    calculation::FixedChargeGroups result;
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
