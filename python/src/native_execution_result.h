#pragma once

#include <chargefw/calculation/calculation.h>

#include <utility>

namespace chargefw::python {

class NativeExecutionResult {
  public:
    explicit NativeExecutionResult(calculation::ExecutionResult result)
        : result_{std::move(result)} {}

    [[nodiscard]] auto result() const noexcept -> const calculation::ExecutionResult& {
        return result_;
    }

  private:
    calculation::ExecutionResult result_;
};

} // namespace chargefw::python
