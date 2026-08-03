#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/core/bond.h>

#include <gemmi/model.hpp>

namespace gemmi::cif {
class Block;
}

#include <vector>

namespace chargefw::adapters::gemmi::bonds {

[[nodiscard]] auto assign(const ::gemmi::Structure& structure, RecordSelection selection,
                          BondStrategy strategy, ::gemmi::cif::Block* mmcif_block = nullptr)
    -> std::vector<core::Bond>;

} // namespace chargefw::adapters::gemmi::bonds
