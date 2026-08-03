#include <chargefw/adapters/gemmi/pdb_input.h>

#include "common_input.h"

#include <gemmi/pdb.hpp>

#include <istream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::adapters::gemmi::pdb_input {

PdbReader::PdbReader(std::istream& input, std::string source,
                     const ::chargefw::adapters::gemmi::RecordSelection selection,
                     const ::chargefw::adapters::gemmi::BondStrategy bond_strategy) {
    std::string contents{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    if (input.bad()) {
        throw std::runtime_error{"failed to read PDB input"};
    }

    const auto structure = ::gemmi::read_pdb_string(contents, source);
    const auto name = structure.name.empty() ? source : structure.name;
    record_ = common_input::make_record(
        structure,
        MoleculeRecordIdentity{.source = std::move(source), .record_index = 0, .record_id = name},
        selection, bond_strategy, nullptr, name);
}

auto PdbReader::next() -> std::optional<ImportedMoleculeRecord> {
    return std::exchange(record_, std::nullopt);
}

} // namespace chargefw::adapters::gemmi::pdb_input
