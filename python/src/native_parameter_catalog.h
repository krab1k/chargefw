#pragma once

#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <utility>
#include <vector>

namespace chargefw::python {

class NativeParameterCatalog {
  public:
    explicit NativeParameterCatalog(std::vector<parameters::ParameterSet> parameter_sets)
        : parameter_sets_{std::move(parameter_sets)} {}

    [[nodiscard]] auto size() const noexcept -> std::size_t {
        return parameter_sets_.size();
    }

    [[nodiscard]] auto parameter_sets() const noexcept
        -> const std::vector<parameters::ParameterSet>& {
        return parameter_sets_;
    }

  private:
    std::vector<parameters::ParameterSet> parameter_sets_;
};

} // namespace chargefw::python
