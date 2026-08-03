#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/molecule_record.h>

#include <istream>
#include <optional>
#include <string>

namespace chargefw::adapters::gemmi::pdb_input {

// Reads a PDB structure through Gemmi. The first PDB MODEL defines atom topology and each model
// becomes a conformer after validating the same atom sequence. Selection can retain all records,
// exclude water, or retain only ATOM records. Alternate locations are excluded: blank locations are
// preferred, otherwise location A, otherwise the first location. BondStrategy::none imports no
// bonds; explicit_bonds imports PDB connectivity; templates adds basic component templates and
// polymer-backbone bonds; hybrid combines explicit_bonds and templates.
class PdbReader {
  public:
    explicit PdbReader(std::istream& input, std::string source = {},
                       ::chargefw::adapters::gemmi::InputOptions options = {});

    [[nodiscard]] auto next() -> std::optional<ImportedMoleculeRecord>;

  private:
    std::optional<ImportedMoleculeRecord> record_;
};

} // namespace chargefw::adapters::gemmi::pdb_input
