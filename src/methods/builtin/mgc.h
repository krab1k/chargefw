#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class MGCMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{.id = "mgc",
                                                 .name = "MGC",
                                                 .full_name = "Molecular Graph Charge",
                                                 .publication = "10.1002/poc.378",
                                                 .priority = 70};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.bond_graph = true;
        requirements.bond_orders = true;
        requirements.element_properties = true;
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