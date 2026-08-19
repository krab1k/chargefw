#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class EQeqCMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "eqeqc",
            .name = "EQeq+C",
            .full_name = "Bond-Order-Corrected Extended Charge Equilibration Method",
            .publication = "10.1021/acs.jctc.5b00037",
            .priority = 140};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.formal_charges = true;
        requirements.element_properties = true;
        requirements.coordinates = true;
        requirements.common_parameters = {"alpha"};
        requirements.atom_parameters = {"Dz"};
        requirements.resources.time = ComplexityTerm::atoms_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_squared;
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
