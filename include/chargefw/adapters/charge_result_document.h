#pragma once

#include <chargefw/adapters/molecule_record.h>
#include <chargefw/charges/charge_collection.h>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace chargefw::adapters {

struct StructuralInputPolicyProvenance {
    std::string selection;
    std::string bonds;
};

// Invocation-wide calculation provenance. The JSON writer serializes the requested inputs and their
// effective resolution as the primary complete result format; other output formats do not consume
// it.
struct RequestedCalculationProvenance {
    std::optional<std::string> method_id;
    std::optional<std::string> parameter_set_id;
    bool permissive_types = false;
    std::optional<std::size_t> full_atom_threshold;
    std::string execution_kind;
    std::optional<double> execution_radius;
    std::optional<std::string> execution_charge_correction;
    std::optional<StructuralInputPolicyProvenance> structural_input_policy;
    std::string conformer_selection;
};

struct EffectiveCalculationProvenance {
    std::optional<std::string> method_id;
    std::optional<std::string> parameter_set_id;
    std::optional<std::string> execution_mode;
    std::optional<double> execution_radius;
    std::optional<std::string> execution_charge_correction;
    std::vector<std::string> warnings;
};

struct CalculationProvenance {
    RequestedCalculationProvenance requested;
    EffectiveCalculationProvenance effective;
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
