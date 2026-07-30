#pragma once

#include <cstdint>

namespace chargefw::adapters::native {

enum class BondFormat : std::uint8_t {
    mol,
    mol2,
};

} // namespace chargefw::adapters::native
