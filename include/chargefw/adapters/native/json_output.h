#pragma once

#include <chargefw/adapters/charge_result_document.h>

#include <iosfwd>

namespace chargefw::adapters::native::json_output {

// Writes a version 1.0 ChargeFW calculation-result JSON document.
class JsonWriter {
  public:
    explicit JsonWriter(std::ostream& output);

    auto write(const ChargeResultDocument& document) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::native::json_output
