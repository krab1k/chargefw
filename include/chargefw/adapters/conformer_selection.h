#pragma once

#include <cstdint>

namespace chargefw::adapters {

enum class ConformerSelection : std::uint8_t {
    first,
    all,
};

} // namespace chargefw::adapters
