#include <chargefw/adapters/gemmi/mmcif_input.h>

#include "common_input.h"

#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>

#include <istream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::adapters::gemmi::mmcif_input {

MmcifReader::MmcifReader(std::istream& input, std::string source,
                         const ::chargefw::adapters::gemmi::RecordSelection selection,
                         const ::chargefw::adapters::gemmi::BondStrategy bond_strategy) {
    std::string contents{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    if (input.bad()) {
        throw std::runtime_error{"failed to read mmCIF input"};
    }

    const auto document =
        ::gemmi::cif::read_memory(contents.data(), contents.size(), source.c_str());
    for (const auto& block : document.blocks) {
        if (!block.has_tag("_atom_site.id")) {
            continue;
        }

        const auto structure = ::gemmi::make_structure_from_block(block);
        records_.push_back(common_input::make_record(
            structure,
            MoleculeRecordIdentity{
                .source = source, .record_index = records_.size(), .record_id = block.name},
            selection, bond_strategy, false, structure.name.empty() ? block.name : structure.name));
    }

    if (records_.empty()) {
        throw std::runtime_error{"mmCIF input contains no coordinate data blocks"};
    }
}

auto MmcifReader::next() -> std::optional<ImportedMoleculeRecord> {
    if (record_index_ == records_.size()) {
        return std::nullopt;
    }

    return std::move(records_[record_index_++]);
}

} // namespace chargefw::adapters::gemmi::mmcif_input
