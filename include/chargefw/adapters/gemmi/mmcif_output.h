#pragma once

#include <chargefw/adapters/charge_result.h>

#include <iosfwd>
#include <string_view>

namespace chargefw::adapters::gemmi::mmcif_output {

// Writes a fresh mmCIF representation of an owned calculation result.
class MmcifWriter {
  public:
    explicit MmcifWriter(std::ostream& output);

    auto write(const ChargeCalculationResult& result, std::string_view generator_name = "ChargeFW",
               std::string_view generator_version = {}) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::gemmi::mmcif_output
