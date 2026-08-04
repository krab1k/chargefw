#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/molecule_record.h>
#include <chargefw/charges/charge_collection.h>

#include <cstdint>
#include <iosfwd>
#include <span>
#include <string_view>

#include <gemmi/cifdoc.hpp>
#include <gemmi/model.hpp>

#include <cstddef>
#include <memory>

namespace chargefw::adapters::gemmi::mmcif_output {

enum class WriteMode : std::uint8_t { replace, append };

struct MmcifSource {
    std::shared_ptr<const ::gemmi::cif::Document> document;
    std::vector<std::size_t> block_indices;
    RecordSelection selection = RecordSelection::all;
};

struct PdbSource {
    ::gemmi::Structure structure;
    RecordSelection selection = RecordSelection::all;
};

// Writes one mmCIF document. Native records are generated as independent UNL data blocks; PDB and
// mmCIF records retain adapter source state so structural content can be converted or enriched.
class MmcifWriter {
  public:
    explicit MmcifWriter(std::ostream& output);

    auto write_generated(std::span<const ImportedMoleculeRecord> records,
                         const charges::ChargeSet& charge_set,
                         std::string_view generator_name = "ChargeFW",
                         std::string_view generator_version = {}) const -> void;

    auto write_pdb(const ImportedMoleculeRecord& record, const charges::ChargeSet& charge_set,
                   const PdbSource& source, std::string_view generator_name = "ChargeFW",
                   std::string_view generator_version = {}) const -> void;

    auto write_mmcif(std::span<const ImportedMoleculeRecord> records,
                     const charges::ChargeSet& charge_set, const MmcifSource& source,
                     std::string_view generator_name = "ChargeFW",
                     std::string_view generator_version = {},
                     WriteMode mode = WriteMode::replace) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::gemmi::mmcif_output
