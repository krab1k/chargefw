#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/features/topology_features.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/methods/method_options.h>

namespace chargefw::methods {
struct CalculationInput {
    const core::Molecule& molecule;
    const features::TopologyFeatures& topology;
    const features::ConformerFeatures* geometry = nullptr;

    const MethodOptions& method_options;
};

} // namespace chargefw::methods
