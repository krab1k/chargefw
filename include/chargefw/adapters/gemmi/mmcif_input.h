#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/molecule_record.h>

#include <gemmi/cifdoc.hpp>

#include <cstddef>
#include <istream>
#include <memory>
#include <optional>
#include <string>

namespace chargefw::adapters::gemmi::mmcif_input {

// Parses an mmCIF document eagerly through Gemmi, then converts one coordinate-bearing data block
// per next() call. Each such block becomes one record; models within a block become conformers
// after atom-sequence validation. Selection, alternate-location, and bond-strategy behavior match
// pdb_input::PdbReader.
class MmcifReader {
  public:
    explicit MmcifReader(std::istream& input, std::string source = {},
                         ::chargefw::adapters::gemmi::InputOptions options = {});

    [[nodiscard]] auto next() -> std::optional<ImportedMoleculeRecord>;
    [[nodiscard]] auto source_document() const -> std::shared_ptr<const ::gemmi::cif::Document>;
    [[nodiscard]] auto source_block_indices() const -> const std::vector<std::size_t>&;
    [[nodiscard]] auto options() const noexcept -> ::chargefw::adapters::gemmi::InputOptions;

  private:
    std::shared_ptr<::gemmi::cif::Document> document_;
    std::string source_;
    ::chargefw::adapters::gemmi::InputOptions options_;
    std::size_t block_index_ = 0;
    std::size_t record_index_ = 0;
    std::vector<std::size_t> source_block_indices_;
};

} // namespace chargefw::adapters::gemmi::mmcif_input
