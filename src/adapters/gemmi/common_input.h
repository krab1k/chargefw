#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/molecule_record.h>

#include <gemmi/model.hpp>

namespace gemmi::cif {
class Block;
}

#include <cstddef>
#include <string>

namespace chargefw::adapters::gemmi::common_input {

[[nodiscard]] auto make_record(const ::gemmi::Structure& structure, MoleculeRecordIdentity identity,
                               RecordSelection selection, BondStrategy bond_strategy,
                               ::gemmi::cif::Block* mmcif_block = nullptr, std::string name = {})
    -> ImportedMoleculeRecord;

} // namespace chargefw::adapters::gemmi::common_input
