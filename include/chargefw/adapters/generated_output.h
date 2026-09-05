#pragma once

#include <chargefw/adapters/molecule_record.h>
#include <chargefw/charges/charge_collection.h>

#include <cstdint>
#include <iosfwd>
#include <span>
#include <string_view>

namespace chargefw::adapters::generated_output {

enum class Format : std::uint8_t { sdf_v2000, sdf_v3000, mol2, mmcif };

// Writes molecular output generated from normalized native records. SDF and MOL2 use the first
// retained conformer; mmCIF writes every retained conformer.
auto write(std::ostream& output, std::span<const ImportedMoleculeRecord> records,
           const charges::ChargeSet& charge_set, Format format,
           std::string_view generator_name = "ChargeFW", std::string_view generator_version = {})
    -> void;

} // namespace chargefw::adapters::generated_output
