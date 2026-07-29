#pragma once

#include <chargefw/adapters/molecule_record.h>

#include <nlohmann/json.hpp>

#include <cstddef>
#include <istream>
#include <optional>
#include <string>

namespace chargefw::adapters::native::json {

// ChargeFW JSON reader for a version 1.0 molecule document. Each call returns the next molecule
// record from the top-level `molecules` array. A null result denotes clean end-of-document; an
// unexpected result denotes one malformed molecule and the next call continues with the next item.
// Molecule array order is source atom/conformer order; atom names are intentionally unsupported.
class JsonReader {
  public:
    explicit JsonReader(std::istream& input, std::string source = {});

    [[nodiscard]] auto next() -> std::optional<::chargefw::adapters::MoleculeRecordResult>;

  private:
    nlohmann::json molecule_records_;
    std::string source_;
    std::size_t record_index_ = 0;
};

} // namespace chargefw::adapters::native::json
