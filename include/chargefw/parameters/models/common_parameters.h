#pragma once

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::parameters {

struct NamedParameter {
    std::string name;
    double value = 0.0;
};

class CommonParameters {
  public:
    CommonParameters() = default;

    explicit CommonParameters(std::vector<NamedParameter> parameters);

    [[nodiscard]] auto parameters() const noexcept -> std::span<const NamedParameter>;

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto empty() const noexcept -> bool;

    [[nodiscard]] auto contains(std::string_view name) const noexcept -> bool;

    [[nodiscard]] auto parameter(std::string_view name) const -> double;

  private:
    std::vector<NamedParameter> parameters_;
};

} // namespace chargefw::parameters