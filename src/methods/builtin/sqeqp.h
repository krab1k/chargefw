#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class SQEqpMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "sqeqp",
            .name = "SQE+qp",
            .full_name = "Split-charge Equilibration with Parameterized Initial Charges",
            .publication = "10.1186/s13321-021-00528-w",
            .priority = 210};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.coordinates = true;
        requirements.atom_parameters = {"electronegativity", "hardness", "width", "q0"};
        requirements.bond_parameters = {"kappa"};
        requirements.resources.time = ComplexityTerm::bonds_cubed;
        requirements.resources.memory = ComplexityTerm::bonds_squared;
        requirements.resources.supports_cutoff = true;
        requirements.resources.supports_cover = true;
        requirements.resources.fragment_target_charge_policy =
            FragmentTargetChargePolicy::proportional_to_atom_count;
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
