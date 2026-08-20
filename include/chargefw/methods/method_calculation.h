#pragma once

#include <chargefw/charges/charge_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>

#include <cstddef>

namespace chargefw::methods {

[[nodiscard]] auto
calculate_charges(const ApplicableMethod& selected,
                  const features::PreparedMoleculeCollection& prepared_collection,
                  std::size_t max_threads = 0) -> charges::ChargeSet;

} // namespace chargefw::methods
