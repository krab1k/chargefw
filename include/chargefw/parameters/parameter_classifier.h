#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/features/topology_features.h>
#include <chargefw/parameters/classification_result.h>
#include <chargefw/parameters/parameter_classification.h>
#include <chargefw/parameters/parameter_set.h>

namespace chargefw::parameters {

[[nodiscard]] auto
try_classify_parameters(const core::Molecule& molecule, const features::TopologyFeatures& topology,
                        const ParameterSet& parameters, const ClassificationOptions& options = {})
    -> ClassificationResult;

[[nodiscard]] auto
classify_parameters(const core::Molecule& molecule, const features::TopologyFeatures& topology,
                    const ParameterSet& parameters, const ClassificationOptions& options = {})
    -> ParameterClassification;

} // namespace chargefw::parameters
