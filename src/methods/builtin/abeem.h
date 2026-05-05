#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class ABEEMMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "abeem",
            .name = "ABEEM",
            .full_name = "Atom-Bond Electronegativity Equalization Method",
            .publication = "10.1021/jp9711048",
            .priority = 190};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.bond_graph = true;
        requirements.formal_charges = true;
        requirements.element_properties = true;
        requirements.coordinates = true;
        requirements.common_parameters = {"k"};
        requirements.atom_parameters = {"a", "b", "c"};
        requirements.bond_parameters = {"A", "B", "C", "D"};
        requirements.resources.time = ComplexityTerm::atoms_plus_bonds_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_plus_bonds_squared;
        requirements.resources.reject_large_without_reduction = true;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const CalculationInput& input) const
        -> charges::AtomicCharges override;
};

} // namespace chargefw::methods::builtin