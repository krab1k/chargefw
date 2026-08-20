#include <chargefw/adapters/gemmi/pdb_input.h>

#include "bonds.h"
#include "selection.h"
#include "structure_import.h"

#include <gemmi/pdb.hpp>

#include <istream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::adapters::gemmi::pdb_input {

PdbReader::PdbReader(std::istream& input, std::string source,
                     const ::chargefw::adapters::gemmi::InputOptions options)
    : options_{options} {
    std::string contents{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    if (input.bad()) {
        throw std::runtime_error{"failed to read PDB input"};
    }

    structure_ = ::gemmi::read_pdb_string(contents, source);
    const auto name = structure_.name.empty() ? source : structure_.name;
    const auto selected = selection::SelectedModel{structure_.models.front(), options.selection};
    auto explicit_bonds = bonds::explicit_pdb(structure_, selected);
    record_ = structure_import::make_record(
        structure_,
        MoleculeRecordIdentity{.source = std::move(source), .record_index = 0, .record_id = name},
        options.selection, options.bond_strategy, options.conformers, std::move(explicit_bonds),
        name);
}

auto PdbReader::next() -> std::optional<ImportedMoleculeRecord> {
    return std::exchange(record_, std::nullopt);
}

auto PdbReader::source_structure() const -> const ::gemmi::Structure& {
    return structure_;
}

auto PdbReader::options() const noexcept -> ::chargefw::adapters::gemmi::InputOptions {
    return options_;
}

} // namespace chargefw::adapters::gemmi::pdb_input
