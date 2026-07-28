#pragma once

#include <chargefw/adapters/molecule_record.h>

#include <expected>
#include <istream>
#include <optional>
#include <string>

namespace chargefw::adapters::native::mol {

// Parses exactly one MOL record. The record ends at M  END; a following SDF separator is not
// consumed. V2000 atom/bond blocks plus M  CHG and the V3000 CTAB atom/bond subset are supported.
[[nodiscard]] auto parse_mol(std::istream& input, MoleculeRecordIdentity identity)
    -> ::chargefw::adapters::MoleculeRecordResult;

// Stateful single-record counterpart to SdfReader. The first call parses the MOL record and later
// calls return clean end-of-file.
class MolReader {
  public:
    explicit MolReader(std::istream& input, std::string source = {});

    [[nodiscard]] auto next() -> std::optional<::chargefw::adapters::MoleculeRecordResult>;

  private:
    std::istream* input_;
    std::string source_;
    bool consumed_ = false;
};

} // namespace chargefw::adapters::native::mol
