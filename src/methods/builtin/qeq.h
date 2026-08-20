#pragma once

#include <chargefw/methods/method.h>

#include <array>

namespace chargefw::methods::builtin {

class QEqMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{.id = "qeq",
                                                 .name = "QEq",
                                                 .full_name = "Charge Equilibration",
                                                 .publication = "10.1021/j100161a070",
                                                 .priority = 170};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.coordinates = true;
        requirements.atom_parameters = {"electronegativity", "hardness"};
        requirements.resources.time = ComplexityTerm::atoms_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_squared;
        requirements.resources.supports_cutoff = true;
        requirements.resources.fragment_target_charge_policy =
            FragmentTargetChargePolicy::proportional_to_atom_count;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const MethodOptionSpec> override {
        static const std::array option_schema{MethodOptionSpec{
            .id = "overlap_term",
            .description = "QEq Coulomb repulsion overlap term",
            .type = MethodOptionType::string,
            .default_value = std::string{"Louwen-Vogt"},
            .choices = {std::string{"Nishimoto-Mataga"}, std::string{"Nishimoto-Mataga-Weiss"},
                        std::string{"Ohno"}, std::string{"Ohno-Klopman"},
                        std::string{"DasGupta-Huzinaga"}, std::string{"Louwen-Vogt"}}}};

        return {option_schema.data(), option_schema.size()};
    }

    [[nodiscard]] auto calculate(const CalculationInput& input) const
        -> charges::AtomicCharges override;
};

} // namespace chargefw::methods::builtin
