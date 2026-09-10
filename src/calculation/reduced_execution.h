#pragma once

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/spatial_fragment.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/classification/parameter_classification.h>

#include <cstddef>
#include <vector>

namespace chargefw::calculation::detail {

struct ReducedChargeContext {
    std::vector<double> reference_charges;
    std::vector<std::vector<std::size_t>> conservation_groups;
};

auto validate_reduced_request(const methods::ApplicableMethod& selected,
                              const ExecutionPolicy& policy, ExecutionMode mode) -> void;
[[nodiscard]] auto prepare_reduced_charge_context(
    const methods::ApplicableMethod& selected, const features::PreparedMolecule& source,
    const parameters::ParameterClassification* source_classification) -> ReducedChargeContext;
[[nodiscard]] auto fragment_target_charge(const ReducedChargeContext& context,
                                          const features::SpatialFragment& fragment) -> double;
auto enforce_conserved_charges(std::vector<double>& values, const ReducedChargeContext& context)
    -> void;
[[nodiscard]] auto
calculate_fragment_charges(const methods::ApplicableMethod& selected,
                           const parameters::ParameterClassification* source_classification,
                           const features::SpatialFragment& fragment,
                           const ReducedChargeContext& charge_context) -> charges::AtomicCharges;

} // namespace chargefw::calculation::detail
