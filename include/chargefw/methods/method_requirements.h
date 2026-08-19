#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace chargefw::methods {

enum class ComplexityTerm : std::uint8_t {
    constant,                 // 1
    atoms,                    // N
    bonds,                    // M
    atoms_plus_bonds,         // N + M
    atoms_squared,            // N^2
    atoms_cubed,              // N^3
    bonds_squared,            // M^2
    bonds_cubed,              // M^3
    atoms_plus_bonds_squared, // (N + M)^2
    atoms_plus_bonds_cubed,   // (N + M)^3
};

enum class FragmentTargetChargePolicy : std::uint8_t {
    unsupported,
    proportional_to_atom_count,
};

struct ResourceRequirements {
    ComplexityTerm time = ComplexityTerm::constant;
    ComplexityTerm memory = ComplexityTerm::constant;

    // Reduced execution is spatial: a method may set either capability only when it also requires
    // coordinates and has a tested executor for the corresponding radius-based approximation.
    bool supports_cutoff = false;
    bool supports_cover = false;
    FragmentTargetChargePolicy fragment_target_charge_policy =
        FragmentTargetChargePolicy::unsupported;
};

struct MethodRequirements {
    bool bond_graph = false;
    bool bond_orders = false;
    bool topological_distances = false;

    bool formal_charges = false;
    bool element_properties = false;

    bool coordinates = false;
    bool spatial_neighbor_search = false;

    std::vector<std::string_view> common_parameters;
    std::vector<std::string_view> atom_parameters;
    std::vector<std::string_view> bond_parameters;

    ResourceRequirements resources{};

    [[nodiscard]] auto requires_common_parameters() const noexcept -> bool {
        return !common_parameters.empty();
    }

    [[nodiscard]] auto requires_atom_parameters() const noexcept -> bool {
        return !atom_parameters.empty();
    }

    [[nodiscard]] auto requires_bond_parameters() const noexcept -> bool {
        return !bond_parameters.empty();
    }

    [[nodiscard]] auto requires_parameters() const noexcept -> bool {
        return requires_common_parameters() || requires_atom_parameters() ||
               requires_bond_parameters();
    }
};

} // namespace chargefw::methods
