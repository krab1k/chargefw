#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class SFKEEMMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "sfkeem",
            .name = "SFKEEM",
            .full_name = "Selfconsistent Functional Kernel Equalized Electronegativity Method",
            .publication = "10.1021/ci050505e",
            .priority = 180};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.formal_charges = true;
        requirements.coordinates = true;
        requirements.common_parameters = {"sigma"};
        requirements.atom_parameters = {"A", "B"};
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
