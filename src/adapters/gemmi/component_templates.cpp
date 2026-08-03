#include "component_templates.h"

#include <algorithm>
#include <iterator>

namespace chargefw::adapters::gemmi::component_templates {
namespace {

#include "component_templates.inc"

struct ComponentTemplateEntry {
    std::string_view component;
    ComponentTemplate component_template;
};

constexpr ComponentTemplateEntry templates[]{
    {"A", {ComponentKind::nucleotide, a_bonds}},
    {"ALA", {ComponentKind::amino_acid, ala_bonds}},
    {"ARG", {ComponentKind::amino_acid, arg_bonds}},
    {"ASN", {ComponentKind::amino_acid, asn_bonds}},
    {"ASP", {ComponentKind::amino_acid, asp_bonds}},
    {"C", {ComponentKind::nucleotide, c_bonds}},
    {"CYS", {ComponentKind::amino_acid, cys_bonds}},
    {"DA", {ComponentKind::nucleotide, da_bonds}},
    {"DC", {ComponentKind::nucleotide, dc_bonds}},
    {"DG", {ComponentKind::nucleotide, dg_bonds}},
    {"DT", {ComponentKind::nucleotide, dt_bonds}},
    {"G", {ComponentKind::nucleotide, g_bonds}},
    {"GLN", {ComponentKind::amino_acid, gln_bonds}},
    {"GLU", {ComponentKind::amino_acid, glu_bonds}},
    {"GLY", {ComponentKind::amino_acid, gly_bonds}},
    {"HIS", {ComponentKind::amino_acid, his_bonds}},
    {"HOH", {ComponentKind::water, hoh_bonds}},
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
    {"U", {ComponentKind::nucleotide, u_bonds}},
    {"VAL", {ComponentKind::amino_acid, val_bonds}},
};

} // namespace

auto find(const std::string_view component) -> std::optional<ComponentTemplate> {
    const auto entry =
        std::lower_bound(std::begin(templates), std::end(templates), component,
                         [](const ComponentTemplateEntry& candidate, const std::string_view id) {
                             return candidate.component < id;
                         });
    if (entry != std::end(templates) && entry->component == component) {
        return entry->component_template;
    }
    return std::nullopt;
}

} // namespace chargefw::adapters::gemmi::component_templates
