#pragma once

#include <chargefw/charges/charge_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <optional>
#include <span>

namespace chargefw::methods {
class Method;
}

namespace chargefw::calculation {

struct CalculationRequest {
    const features::PreparedMoleculeCollection& molecules;
    std::span<const methods::Method* const> candidate_methods;
    std::span<const parameters::ParameterSet> parameter_sets;
};

struct CalculationResult {
    std::optional<charges::ChargeSet> charges;
    methods::ApplicabilityResult applicability;

    [[nodiscard]] auto calculated() const noexcept -> bool {
        return charges.has_value();
    }
};

// Selects the applicable candidate with the highest method priority, then the highest parameter-set
// priority. Equal priorities are resolved by method ID, then parameter-set ID, in lexicographic
// order.
[[nodiscard]] auto calculate(const CalculationRequest& request) -> CalculationResult;

} // namespace chargefw::calculation
