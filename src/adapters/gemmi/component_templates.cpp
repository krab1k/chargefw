#include "component_templates.h"

#include <algorithm>
#include <array>

namespace chargefw::adapters::gemmi::component_templates {
namespace {

#include "component_templates.inc"

struct ComponentTemplateEntry {
    std::string_view component;
    ComponentTemplate component_template;
};

constexpr std::array<ComponentTemplateEntry, 20> amino_acid_templates{{
    {"ALA", {ComponentKind::amino_acid, ala_bonds}},
    {"ARG", {ComponentKind::amino_acid, arg_bonds}},
    {"ASN", {ComponentKind::amino_acid, asn_bonds}},
    {"ASP", {ComponentKind::amino_acid, asp_bonds}},
    {"CYS", {ComponentKind::amino_acid, cys_bonds}},
    {"GLN", {ComponentKind::amino_acid, gln_bonds}},
    {"GLU", {ComponentKind::amino_acid, glu_bonds}},
    {"GLY", {ComponentKind::amino_acid, gly_bonds}},
    {"HIS", {ComponentKind::amino_acid, his_bonds}},
    {"ILE", {ComponentKind::amino_acid, ile_bonds}},
    {"LEU", {ComponentKind::amino_acid, leu_bonds}},
    {"LYS", {ComponentKind::amino_acid, lys_bonds}},
    {"MET", {ComponentKind::amino_acid, met_bonds}},
    {"PHE", {ComponentKind::amino_acid, phe_bonds}},
    {"PRO", {ComponentKind::amino_acid, pro_bonds}},
    {"SER", {ComponentKind::amino_acid, ser_bonds}},
    {"THR", {ComponentKind::amino_acid, thr_bonds}},
    {"TRP", {ComponentKind::amino_acid, trp_bonds}},
    {"TYR", {ComponentKind::amino_acid, tyr_bonds}},
    {"VAL", {ComponentKind::amino_acid, val_bonds}},
}};

constexpr std::array<ComponentTemplateEntry, 8> nucleotide_templates{{
    {"A", {ComponentKind::nucleotide, a_bonds}},
    {"C", {ComponentKind::nucleotide, c_bonds}},
    {"DA", {ComponentKind::nucleotide, da_bonds}},
    {"DC", {ComponentKind::nucleotide, dc_bonds}},
    {"DG", {ComponentKind::nucleotide, dg_bonds}},
    {"DT", {ComponentKind::nucleotide, dt_bonds}},
    {"G", {ComponentKind::nucleotide, g_bonds}},
    {"U", {ComponentKind::nucleotide, u_bonds}},
}};

constexpr std::array<ComponentTemplateEntry, 1> water_templates{{
    {"HOH", {ComponentKind::water, hoh_bonds}},
}};

template <std::size_t Size>
[[nodiscard]] auto find_in(const std::array<ComponentTemplateEntry, Size>& templates,
                           const std::string_view component) -> std::optional<ComponentTemplate> {
    const auto entry =
        std::lower_bound(templates.begin(), templates.end(), component,
                         [](const ComponentTemplateEntry& candidate, const std::string_view id) {
                             return candidate.component < id;
                         });
    if (entry != templates.end() && entry->component == component) {
        return entry->component_template;
    }
    return std::nullopt;
}

} // namespace

auto find(const std::string_view component) -> std::optional<ComponentTemplate> {
    if (const auto component_template = find_in(amino_acid_templates, component)) {
        return component_template;
    }
    if (const auto component_template = find_in(nucleotide_templates, component)) {
        return component_template;
    }
    return find_in(water_templates, component);
}

} // namespace chargefw::adapters::gemmi::component_templates
