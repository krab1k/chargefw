#pragma once

#include <chargefw/adapters/molecule_record.h>
#include <chargefw/charges/charge_collection.h>

#include <optional>
#include <string>
#include <vector>

namespace chargefw::adapters {

struct ChargeResultRecord {
    MoleculeRecordIdentity identity;
    MoleculeRecordMapping mapping;
    std::optional<charges::ChargeSet> charges;
};

struct ChargeResultDocument {
    std::string generator_name;
    std::string generator_version;
    std::vector<ChargeResultRecord> records;
};

} // namespace chargefw::adapters
