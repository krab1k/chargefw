#pragma once

#include <chargefw/methods/method.h>

#include <span>
#include <vector>

namespace chargefw::methods::builtin {

namespace sqe_core {

[[nodiscard]] auto calculate(const CalculationInput& input,
                             std::span<const double> initial_charges = {}) -> std::vector<double>;

} // namespace sqe_core

class SQEMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{.id = "sqe",
                                                 .name = "SQE",
                                                 .full_name = "Split-charge Equilibration",
                                                 .publication = "10.1063/1.2346671",
                                                 .priority = 90};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.coordinates = true;
        requirements.atom_parameters = {"electronegativity", "hardness", "width"};
        requirements.bond_parameters = {"kappa"};
        requirements.resources.time = ComplexityTerm::bonds_cubed;
        requirements.resources.memory = ComplexityTerm::bonds_squared;
        requirements.resources.supports_cutoff = true;
        requirements.resources.supports_cover = true;
        requirements.resources.fragment_target_charge_policy = FragmentTargetChargePolicy::zero;
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
