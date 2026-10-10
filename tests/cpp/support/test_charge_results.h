#pragma once

#include <chargefw/adapters/charge_result.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/charge_collection.h>

#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::test {

/// Wraps a charge set as a full-execution result owned by the given records.
[[nodiscard]] inline auto calculation_result(std::vector<adapters::ImportedMoleculeRecord> records,
                                             charges::ChargeSet charge_set)
    -> adapters::ChargeCalculationResult {
    const auto method_id = std::string{charge_set.method_id()};
    const auto parameter_set_id = charge_set.parameter_set_id().transform(
        [](const std::string_view value) { return std::string{value}; });
    return adapters::make_charge_calculation_result(
        std::move(records), {},
        {.charges = std::move(charge_set),
         .effective = calculation::EffectiveCalculation{.method_id = method_id,
                                                        .parameter_set_id = parameter_set_id,
                                                        .execution_policy =
                                                            calculation::ExecutionPolicy{}}});
}

} // namespace chargefw::test
