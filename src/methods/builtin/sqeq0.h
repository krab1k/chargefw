#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class SQEq0Method final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "sqeq0",
            .name = "SQE+q0",
            .full_name = "Split-charge Equilibration with Initial Formal Charges",
            .publication = "10.1021/ct200512e",
            .priority = 80};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.bond_graph = true;
        requirements.formal_charges = true;
        requirements.coordinates = true;
        requirements.atom_parameters = {"electronegativity", "hardness", "width"};
        requirements.bond_parameters = {"kappa"};
        requirements.resources.time = ComplexityTerm::bonds_cubed;
        requirements.resources.memory = ComplexityTerm::bonds_squared;
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
