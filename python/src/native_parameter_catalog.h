#pragma once

#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <string>
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

    [[nodiscard]] auto descriptors() const -> std::vector<parameters::ParameterSetMetadata> {
        std::vector<parameters::ParameterSetMetadata> result;
        result.reserve(parameter_sets_.size());
        for (const auto& parameter_set : parameter_sets_) {
            result.push_back(parameter_set.metadata());
        }
        return result;
    }

  private:
    std::vector<parameters::ParameterSet> parameter_sets_;
};

} // namespace chargefw::python
