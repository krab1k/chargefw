#pragma once

#include <chargefw/adapters/molecule_record.h>

#include <cstddef>
#include <istream>
#include <optional>
#include <string>

namespace chargefw::adapters::native::sdf_input {

// Bounded-memory SDF reader. Each call consumes at most one SDF record. A null result denotes clean
// end-of-file; an unexpected result denotes one malformed record and the next call starts later.
class SdfReader {
  public:
    explicit SdfReader(std::istream& input, std::string source = {});

    [[nodiscard]] auto next() -> std::optional<::chargefw::adapters::MoleculeRecordResult>;

  private:
    std::istream* input_;
    std::string source_;
    std::size_t record_index_ = 0;
};

} // namespace chargefw::adapters::native::sdf_input
