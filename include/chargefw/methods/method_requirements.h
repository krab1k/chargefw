#pragma once

#include <cstddef>
#include <string_view>
#include <vector>

namespace chargefw::methods {

inline constexpr std::size_t default_large_molecule_atom_threshold = 20'000;

enum class ComplexityTerm {
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

struct ResourceRequirements {
    ComplexityTerm time = ComplexityTerm::constant;
    ComplexityTerm memory = ComplexityTerm::constant;

    bool supports_cutoff = false;
    bool supports_cover = false;

    std::size_t large_molecule_atom_threshold = default_large_molecule_atom_threshold;
    bool reject_large_without_reduction = false;
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