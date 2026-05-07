#pragma once

#include <chargefw/methods/method.h>

#include <array>

namespace chargefw::methods::builtin {

class DENRMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{.id = "denr",
                                                 .name = "DENR",
                                                 .full_name =
                                                     "Dynamical Electronegativity Relaxation",
                                                 .publication = "10.1080/10629360701844142",
                                                 .priority = 50};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.bond_graph = true;
        requirements.atom_parameters = {"electronegativity", "hardness"};
        requirements.resources.time = ComplexityTerm::atoms_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_squared;
        requirements.resources.reject_large_without_reduction = true;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const MethodOptionSpec> override {
        static const std::array option_schema{
            MethodOptionSpec{.id = "step",
                             .description = "DENR relaxation step size",
                             .type = MethodOptionType::floating_point,
                             .default_value = 0.1,
                             .choices = {}},
            MethodOptionSpec{.id = "iterations",
                             .description = "Number of DENR iterations",
                             .type = MethodOptionType::integer,
                             .default_value = 3,
                             .choices = {}}};

        return {option_schema.data(), option_schema.size()};
    }

    [[nodiscard]] auto calculate(const CalculationInput& input) const
        -> charges::AtomicCharges override;
};

} // namespace chargefw::methods::builtin