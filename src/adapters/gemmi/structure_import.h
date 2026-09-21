#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/molecule_record.h>

#include <gemmi/model.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace chargefw::adapters::gemmi::structure_import {

struct SourceModelMapping {
    SourceConformerReference conformer;
};

[[nodiscard]] auto make_record(const ::gemmi::Structure& structure, MoleculeRecordIdentity identity,
                               RecordSelection selection, BondStrategy bond_strategy,
                               ConformerSelection conformer_selection = ConformerSelection::all,
                               std::vector<core::Bond> explicit_bonds = {}, std::string name = {},
                               MolecularSourceFormat format = MolecularSourceFormat::pdb,
                               std::vector<SourceModelMapping> source_models = {},
                               SourceConnectivity source_connectivity = SourceConnectivity::absent)
    -> ImportedMoleculeRecord;

} // namespace chargefw::adapters::gemmi::structure_import
