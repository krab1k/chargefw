#pragma once

#include <chargefw/methods/method.h>

#include <array>

namespace chargefw::methods::builtin {

class MPEOEMethod final : public Method {
public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "mpeoe",
            .name = "MPEOE",
            .full_name = "Modified Partial Equalization of Atomic Electronegativity",
            .publication = "10.1021/j100374a066",
            .priority = 110
        };

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.bond_graph = true;
        requirements.common_parameters = {"Hplus"};
        requirements.atom_parameters = {"A", "B"};
        requirements.bond_parameters = {"f"};
        requirements.resources.time = ComplexityTerm::atoms_plus_bonds;
        requirements.resources.memory = ComplexityTerm::atoms;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const MethodOptionSpec> override {
        return {
            option_schema_.data(),
            option_schema_.size()
        };
    }

    [[nodiscard]] auto calculate(const CalculationInput& input) const
        -> charges::AtomicCharges override;

private:
    inline static const std::array<MethodOptionSpec, 1> option_schema_{
        MethodOptionSpec{
            .id = "iters",
            .description = "Number of MPEOE iterations",
            .type = MethodOptionType::integer,
            .default_value = 7,
            .choices = {}
        }
    };
};

} // namespace chargefw::methods::builtin