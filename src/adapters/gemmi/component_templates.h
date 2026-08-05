#pragma once

#include <chargefw/core/bond.h>

#include <cstdint>
#include <span>
#include <string_view>

namespace chargefw::adapters::gemmi::component_templates {

struct BondTemplate {
    std::string_view first;
    std::string_view second;
    core::BondOrder order;
};

enum class ComponentKind : std::uint8_t {
    amino_acid,
    nucleotide,
    water,
};

struct ComponentTemplate {
    ComponentKind kind;
    std::span<const BondTemplate> bonds;
};

[[nodiscard]] auto find(std::string_view component) -> const ComponentTemplate*;

} // namespace chargefw::adapters::gemmi::component_templates
