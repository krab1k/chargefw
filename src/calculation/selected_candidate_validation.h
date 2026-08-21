#pragma once

#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>

namespace chargefw::calculation::detail {

auto validate_selected_candidate(const methods::ApplicableMethod& selected,
                                 const features::PreparedMoleculeCollection& molecules) -> void;

} // namespace chargefw::calculation::detail
