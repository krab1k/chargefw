#pragma once

#include <chargefw/methods/method.h>

#include <string_view>
#include <vector>

namespace chargefw::methods::builtin {

namespace eqeq_core {

[[nodiscard]] auto calculate(const CalculationInput& input, std::string_view method_name)
    -> std::vector<double>;

} // namespace eqeq_core

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
        requirements.coordinates = true;
        requirements.resources.time = ComplexityTerm::atoms_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_squared;
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
