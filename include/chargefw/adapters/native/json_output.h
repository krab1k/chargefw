#pragma once

#include <chargefw/adapters/charge_result_document.h>

#include <iosfwd>
#include <optional>
#include <string_view>

namespace chargefw::adapters::native::json_output {

// Writes a version 1.0 ChargeFW calculation-result JSON document.
class JsonWriter {
  public:
    explicit JsonWriter(std::ostream& output);

    auto write(const ChargeResultDocument& document) const -> void;
    auto write(const ChargeCalculationResult& result, std::string_view generator_name,
               std::string_view generator_version,
               std::optional<ExecutionMetrics> execution_metrics = std::nullopt) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::native::json_output
