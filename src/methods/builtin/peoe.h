#pragma once

#include <chargefw/methods/method.h>

#include <array>

namespace chargefw::methods::builtin {

class PEOEMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "peoe",
            .name = "PEOE",
            .full_name = "Partial Equalization of Atomic Electronegativity",
            .publication = "10.1016/0040-4020(80)80168-2",
            .priority = 120};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.common_parameters = {"dampH"};
        requirements.atom_parameters = {"A", "B", "C"};
        requirements.resources.time = ComplexityTerm::atoms_plus_bonds;
        requirements.resources.memory = ComplexityTerm::atoms;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const MethodOptionSpec> override {
        static const std::array option_schema{
            MethodOptionSpec{.id = "iters",
                             .description = "Number of PEOE iterations",
                             .type = MethodOptionType::integer,
                             .default_value = 7,
                             .choices = {}}};
        return {option_schema.data(), option_schema.size()};
    }

    [[nodiscard]] auto calculate(const CalculationInput& input) const
        -> charges::AtomicCharges override;
};

} // namespace chargefw::methods::builtin
