#pragma once

#include <chargefw/methods/method_options.h>
#include <chargefw/methods/feature_requirements.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/charges/atomic_charges.h>

#include <span>
#include <string_view>

namespace chargefw::methods {

class Method {
  public:
    virtual ~Method() = default;

    [[nodiscard]] virtual auto metadata() const noexcept -> const MethodMetadata& = 0;

    [[nodiscard]] auto id() const noexcept -> std::string_view {
        return metadata().id;
    }

    [[nodiscard]] virtual auto requirements() const noexcept -> FeatureRequirements = 0;

    [[nodiscard]] virtual auto option_schema() const noexcept
        -> std::span<const MethodOptionSpec> = 0;

    [[nodiscard]] virtual auto calculate(const CalculationInput& input) const
        -> charges::AtomicCharges = 0;
};

} // namespace chargefw::methods
