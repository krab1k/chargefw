#pragma once

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/molecule.h>
#include <chargefw/features/spatial_fragment.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/classification/parameter_classification.h>

#include <optional>
#include <string>
#include <vector>

namespace chargefw::calculation::detail {

auto validate_reduced_request(const methods::ApplicableMethod& selected,
                              const ExecutionPolicy& policy, ExecutionMode mode) -> void;
[[nodiscard]] auto fragment_target_charge(methods::FragmentTargetChargePolicy policy,
                                          const core::Molecule& source,
                                          const features::SpatialFragment& fragment) -> double;
[[nodiscard]] auto final_target_charge(methods::FragmentTargetChargePolicy policy,
                                       const core::Molecule& source) -> double;
auto validate_fragment_charges(const methods::Method& method,
                               const features::SpatialFragment& fragment,
                               const charges::AtomicCharges& charges) -> void;
auto apply_charge_correction(std::vector<double>& values, double target_charge,
                             ChargeCorrectionPolicy policy) -> void;
[[nodiscard]] auto parameter_set_id_for(const methods::ApplicableMethod& selected)
    -> std::optional<std::string>;
[[nodiscard]] auto
calculate_fragment_charges(const methods::ApplicableMethod& selected, const core::Molecule& source,
                           const parameters::ParameterClassification* source_classification,
                           const features::SpatialFragment& fragment) -> charges::AtomicCharges;

} // namespace chargefw::calculation::detail
