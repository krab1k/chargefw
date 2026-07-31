#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/molecule_record.h>

#include <gemmi/model.hpp>

#include <cstddef>
#include <string>

namespace chargefw::adapters::gemmi::common_input {

[[nodiscard]] auto make_record(const ::gemmi::Structure& structure, MoleculeRecordIdentity identity,
                               RecordSelection selection, std::string name = {})
    -> ImportedMoleculeRecord;

} // namespace chargefw::adapters::gemmi::common_input
