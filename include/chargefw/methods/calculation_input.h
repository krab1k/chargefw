#pragma once

#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_options.h>

namespace chargefw::methods {

struct CalculationInput {
    const features::PreparedMolecule& prepared_molecule;
    const features::ConformerFeatures* geometry = nullptr;

    const MethodOptions& method_options;
};

} // namespace chargefw::methods