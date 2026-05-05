#include <chargefw/core/position.h>

#include <cmath>

namespace chargefw::core {

auto squared_distance(const Position& first, const Position& second) noexcept -> double {
    const auto dx = first.x - second.x;
    const auto dy = first.y - second.y;
    const auto dz = first.z - second.z;

    return dx * dx + dy * dy + dz * dz;
}

auto distance(const Position& first, const Position& second) noexcept -> double {
    return std::sqrt(squared_distance(first, second));
}

} // namespace chargefw::core