#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/molecule_record.h>

#include <gemmi/model.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace chargefw::adapters::gemmi::structure_import {

[[nodiscard]] auto make_record(const ::gemmi::Structure& structure, MoleculeRecordIdentity identity,
                               RecordSelection selection, BondStrategy bond_strategy,
                               std::vector<core::Bond> explicit_bonds = {}, std::string name = {})
    -> ImportedMoleculeRecord;

} // namespace chargefw::adapters::gemmi::structure_import
