#pragma once

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>

namespace chargefw::calculation {

[[nodiscard]] auto calculate_cutoff_charges(const methods::ApplicableMethod& selected,
                                            const features::PreparedMoleculeCollection& molecules,
                                            const ExecutionPolicy& policy,
                                            std::size_t max_threads = 0) -> charges::ChargeSet;

} // namespace chargefw::calculation
