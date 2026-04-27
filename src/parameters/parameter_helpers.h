#pragma once

#include <chargefw/parameters/common_parameters.h>

#include <span>
#include <string_view>

namespace chargefw::parameters::detail {

auto validate_named_parameters(std::span<const NamedParameter> parameters, std::string_view context)
    -> void;

[[nodiscard]] auto contains_named_parameter(std::span<const NamedParameter> parameters,
                                            std::string_view name) noexcept -> bool;

[[nodiscard]] auto named_parameter(std::span<const NamedParameter> parameters,
                                   std::string_view name, std::string_view context) -> double;

} // namespace chargefw::parameters::detail