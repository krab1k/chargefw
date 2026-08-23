#pragma once

#include <chargefw/calculation/observer.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>

#include <cstddef>

namespace chargefw::calculation {

[[nodiscard]] auto calculate_full_charges(const methods::ApplicableMethod& selected,
                                          const features::PreparedMoleculeCollection& molecules,
                                          std::size_t max_threads,
                                          const CalculationObserver& observer)
    -> charges::ChargeSet;

} // namespace chargefw::calculation
