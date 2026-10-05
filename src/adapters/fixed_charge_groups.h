#pragma once

#include <chargefw/adapters/molecule_record.h>
#include <chargefw/calculation/assessment.h>

#include <span>
#include <string>
#include <vector>

namespace chargefw::adapters::detail {

[[nodiscard]] auto fixed_charge_ion_names(bool common_only) -> std::vector<std::string>;

[[nodiscard]] auto resolve_fixed_charge_groups(std::span<const ImportedMoleculeRecord> records,
                                               std::span<const std::string> component_ids)
    -> calculation::FixedChargeGroups;

} // namespace chargefw::adapters::detail
