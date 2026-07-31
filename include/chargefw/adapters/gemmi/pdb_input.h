#pragma once

#include <chargefw/adapters/molecule_record.h>

#include <istream>
#include <optional>
#include <string>

namespace chargefw::adapters::gemmi::pdb_input {

enum class RecordSelection {
    all,
    polymers_and_ligands,
    polymers,
};

// Reads a PDB structure through Gemmi. The first PDB MODEL defines atom topology and each model
// becomes a conformer after validating the same atom sequence. Selection can retain all records,
// exclude water, or retain only ATOM records. Alternate locations are excluded: blank locations are
// preferred, otherwise location A, otherwise the first location. Bonds are omitted because PDB
// connectivity requires an explicit, future bond-adding policy.
class PdbReader {
  public:
    explicit PdbReader(std::istream& input, std::string source = {},
                       RecordSelection selection = RecordSelection::all);

    [[nodiscard]] auto next() -> std::optional<ImportedMoleculeRecord>;

  private:
    std::optional<ImportedMoleculeRecord> record_;
};

} // namespace chargefw::adapters::gemmi::pdb_input
