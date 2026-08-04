#pragma once

#include <chargefw/adapters/molecule_record.h>

#include <cstddef>
#include <istream>
#include <optional>
#include <string>

namespace chargefw::adapters::native::mol2_input {

// Bounded-memory Tripos MOL2 reader. Each call consumes one @<TRIPOS>MOLECULE record. The reader
// supports MOLECULE, ATOM, and BOND sections with standard atom types and numeric bonds. Aromatic
// bond types are imported as single bonds.
class Mol2Reader {
  public:
    explicit Mol2Reader(std::istream& input, std::string source = {});

    [[nodiscard]] auto next() -> std::optional<ImportedMoleculeRecord>;

  private:
    std::istream* input_;
    std::string source_;
    std::size_t record_index_ = 0;
};

} // namespace chargefw::adapters::native::mol2_input
