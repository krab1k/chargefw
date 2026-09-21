#pragma once

#include <chargefw/adapters/charge_result.h>

#include <utility>

namespace chargefw::python {

class NativeExecutionResult {
  public:
    explicit NativeExecutionResult(adapters::ChargeCalculationResult result)
        : result_{std::move(result)} {}

    [[nodiscard]] auto result() const noexcept -> const adapters::ChargeCalculationResult& {
        return result_;
    }

  private:
    adapters::ChargeCalculationResult result_;
};

} // namespace chargefw::python
