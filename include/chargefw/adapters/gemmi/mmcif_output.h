#pragma once

#include <chargefw/adapters/charge_result.h>

#include <iosfwd>
#include <string_view>

#include <gemmi/cifdoc.hpp>

namespace chargefw::adapters::gemmi::mmcif_output {

// Writes fresh result documents or strictly annotates an unchanged caller-owned source document.
class MmcifWriter {
  public:
    explicit MmcifWriter(std::ostream& output);

    auto write(const ChargeCalculationResult& result, std::string_view generator_name = "ChargeFW",
               std::string_view generator_version = {}) const -> void;

    auto write_attached(const ChargeCalculationResult& result, const ::gemmi::cif::Document& source,
                        bool overwrite = false, std::string_view generator_name = "ChargeFW",
                        std::string_view generator_version = {}) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::gemmi::mmcif_output
