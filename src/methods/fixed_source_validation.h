#pragma once

#include <chargefw/methods/calculation_input.h>

#include <string_view>

namespace chargefw::methods::detail {

auto validate_fixed_source_input(const CalculationInput& input, std::string_view method_name)
    -> void;

} // namespace chargefw::methods::detail
