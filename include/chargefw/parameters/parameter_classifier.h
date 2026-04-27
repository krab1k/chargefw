#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/features/topology_features.h>
#include <chargefw/parameters/parameter_classification.h>
#include <chargefw/parameters/parameter_set.h>

namespace chargefw::parameters {

[[nodiscard]] auto classify_parameters(const core::Molecule& molecule,
                                       const features::TopologyFeatures& topology,
                                       const ParameterSet& parameters) -> ParameterClassification;

} // namespace chargefw::parameters