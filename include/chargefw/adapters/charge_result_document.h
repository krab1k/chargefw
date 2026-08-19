#pragma once

#include <chargefw/adapters/molecule_record.h>
#include <chargefw/charges/charge_collection.h>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace chargefw::adapters {

// Invocation-wide calculation provenance. The JSON writer serializes this as the primary complete
// result format; other output formats intentionally do not consume it yet.
struct CalculationProvenance {
    std::optional<std::string> effective_execution_mode;
    std::optional<double> radius;
    std::optional<std::string> charge_correction;
    bool permissive_types = false;
    std::optional<std::size_t> full_atom_threshold;
    std::vector<std::string> warnings;
};

struct ChargeResultRecord {
    MoleculeRecordIdentity identity;
    MoleculeRecordMapping mapping;
    std::optional<charges::ChargeSet> charges;
};

struct ChargeResultDocument {
    std::string generator_name;
    std::string generator_version;
    std::vector<ChargeResultRecord> records;
    std::optional<CalculationProvenance> calculation_provenance;
};

} // namespace chargefw::adapters
