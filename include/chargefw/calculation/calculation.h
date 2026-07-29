#pragma once

#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <optional>
#include <span>
#include <string>
#include <vector>

namespace chargefw::methods {
class Method;
}

namespace chargefw::calculation {

struct CalculationRequest {
    // Non-owning input; the referenced collection must outlive calculate().
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
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

// Owns application-facing inputs so language bindings and adapters do not need to manage native
// method pointers, parameter spans, or prepared-feature lifetimes. Method-specific options remain
// available through the low-level CalculationRequest for advanced native callers.
struct ApplicationCalculationRequest {
    core::MoleculeCollection molecules;
    std::vector<parameters::ParameterSet> parameter_sets;
    std::optional<std::string> method_id;
    std::optional<std::string> parameter_set_id;
};

using ApplicationCalculationResult = CalculationResult;

// Selects the applicable candidate with the highest method priority, then the highest parameter-set
// priority. Equal priorities are resolved by method ID, then parameter-set ID, in lexicographic
// order.
[[nodiscard]] auto calculate(const CalculationRequest& request) -> CalculationResult;

// Calculates an owned molecule collection using registered methods and caller-owned parameter data.
// Omitted IDs enable deterministic automatic selection. Specified IDs restrict selection to the
// exact method or parameter set; unavailable or inapplicable selections are reported as errors.
[[nodiscard]] auto calculate(const ApplicationCalculationRequest& request)
    -> ApplicationCalculationResult;

} // namespace chargefw::calculation
