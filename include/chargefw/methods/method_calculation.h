#pragma once

#include <chargefw/charges/charge_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>

namespace chargefw::methods {

[[nodiscard]] auto calculate_charges(const ApplicableMethod& selected,
                                     const features::PreparedMoleculeCollection& molecules)
    -> charges::ChargeSet;

} // namespace chargefw::methods