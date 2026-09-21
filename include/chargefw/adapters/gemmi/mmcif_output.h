#pragma once

#include <chargefw/adapters/charge_result.h>
#include <chargefw/adapters/molecule_record.h>
#include <chargefw/charges/charge_collection.h>

#include <iosfwd>
#include <span>
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

    auto write_generated(std::span<const ImportedMoleculeRecord> records,
                         const charges::ChargeSet& charge_set,
                         std::string_view generator_name = "ChargeFW",
                         std::string_view generator_version = {}) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::gemmi::mmcif_output
