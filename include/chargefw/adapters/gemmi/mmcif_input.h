#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/molecule_record.h>

#include <istream>
#include <optional>
#include <string>
#include <vector>

namespace chargefw::adapters::gemmi::mmcif_input {

// Reads coordinate-bearing mmCIF data blocks through Gemmi. Each data block becomes one record;
// models within a block become conformers after atom-sequence validation. Selection, alternate
// location, and deferred-bond behavior match pdb_input::PdbReader.
class MmcifReader {
  public:
    explicit MmcifReader(std::istream& input, std::string source = {},
                         ::chargefw::adapters::gemmi::RecordSelection selection =
                             ::chargefw::adapters::gemmi::RecordSelection::all);

    [[nodiscard]] auto next() -> std::optional<ImportedMoleculeRecord>;

  private:
    std::vector<ImportedMoleculeRecord> records_;
    std::size_t record_index_ = 0;
};

} // namespace chargefw::adapters::gemmi::mmcif_input
