#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class EQeqMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{.id = "eqeq",
                                                 .name = "EQeq",
                                                 .full_name =
                                                     "Extended Charge Equilibration Method",
                                                 .publication = "10.1021/jz3008485",
                                                 .priority = 150};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.formal_charges = true;
        requirements.element_properties = true;
        requirements.coordinates = true;
        requirements.resources.time = ComplexityTerm::atoms_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_squared;
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