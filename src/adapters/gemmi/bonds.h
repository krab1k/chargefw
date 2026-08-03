#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/core/bond.h>

#include "selection.h"

#include <gemmi/model.hpp>

namespace gemmi::cif {
class Block;
}

#include <vector>

namespace chargefw::adapters::gemmi::bonds {

[[nodiscard]] auto assign(const selection::SelectedModel& model, BondStrategy strategy,
                          std::vector<core::Bond> explicit_bonds = {}) -> std::vector<core::Bond>;

[[nodiscard]] auto explicit_pdb(const ::gemmi::Structure& structure,
                                const selection::SelectedModel& model) -> std::vector<core::Bond>;

[[nodiscard]] auto explicit_mmcif(const ::gemmi::Structure& structure, ::gemmi::cif::Block& block,
                                  const selection::SelectedModel& model) -> std::vector<core::Bond>;

} // namespace chargefw::adapters::gemmi::bonds
