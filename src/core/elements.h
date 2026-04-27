#pragma once

#include <chargefw/core/periodic_table.h>

#include <span>

namespace chargefw::core::detail {

[[nodiscard]] auto elements() noexcept -> std::span<const Element>;

} // namespace chargefw::core::detail
