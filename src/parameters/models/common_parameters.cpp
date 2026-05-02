#include <chargefw/parameters/models/common_parameters.h>

#include "parameters/parameter_helpers.h"

#include <utility>

namespace chargefw::parameters {

CommonParameters::CommonParameters(std::vector<NamedParameter> parameters)
    : parameters_{std::move(parameters)} {
    detail::validate_named_parameters(parameters_, "common");
}

auto CommonParameters::parameters() const noexcept -> std::span<const NamedParameter> {
    return parameters_;
}

auto CommonParameters::size() const noexcept -> std::size_t {
    return parameters_.size();
}

auto CommonParameters::empty() const noexcept -> bool {
    return parameters_.empty();
}

auto CommonParameters::contains(const std::string_view name) const noexcept -> bool {
    return detail::contains_named_parameter(parameters_, name);
}

auto CommonParameters::parameter(const std::string_view name) const -> double {
    return detail::named_parameter(parameters_, name, "common");
}

} // namespace chargefw::parameters