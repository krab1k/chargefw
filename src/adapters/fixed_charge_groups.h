#pragma once

#include <chargefw/adapters/molecule_record.h>
#include <chargefw/calculation/assessment.h>

#include <span>
#include <string>

namespace chargefw::adapters::detail {

[[nodiscard]] auto resolve_fixed_charge_groups(std::span<const ImportedMoleculeRecord> records,
                                               std::span<const std::string> component_ids)
    -> calculation::FixedChargeEmbedding;

} // namespace chargefw::adapters::detail
