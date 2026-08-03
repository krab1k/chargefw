#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/core/bond.h>

#include <gemmi/model.hpp>

#include <vector>

namespace chargefw::adapters::gemmi::bonds {

[[nodiscard]] auto assign_template_bonds(const ::gemmi::Model& model, RecordSelection selection)
    -> std::vector<core::Bond>;

} // namespace chargefw::adapters::gemmi::bonds
