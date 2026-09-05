#pragma once

#include <chargefw/adapters/molecule_record.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/methods/method_options.h>

#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::adapters {

enum class DiagnosticSeverity : std::uint8_t { info, warning, error };

struct ResultDiagnostic {
    DiagnosticSeverity severity = DiagnosticSeverity::error;
    std::string code;
    std::string message;
    // Structured indices are zero-based. Source line numbers are one-based.
    std::optional<std::size_t> molecule_index;
    std::optional<std::size_t> atom_index;
    std::optional<std::size_t> bond_index;
    std::optional<std::size_t> conformer_index;
    std::optional<std::size_t> line;
};

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
    std::optional<std::size_t> cutoff_atom_threshold;
    std::optional<std::size_t> cover_atom_threshold;
    std::size_t max_threads = 0;
    std::string execution_kind;
    std::optional<double> execution_radius;
    std::optional<std::string> execution_charge_correction;
    std::optional<StructuralInputPolicyProvenance> structural_input_policy;
    std::optional<std::string> conformer_selection;
    std::map<std::string, methods::MethodOptions> method_options;
};

struct EffectiveCalculationProvenance {
    std::optional<std::string> method_id;
    std::optional<std::string> parameter_set_id;
    std::optional<std::string> execution_mode;
    std::optional<double> execution_radius;
    std::optional<std::string> execution_charge_correction;
    std::vector<std::string> warnings;
    std::map<std::string, methods::MethodOptions> method_options;
};

struct ExecutionMetrics {
    std::string started_at;
    std::string ended_at;
    double runtime_seconds = 0.0;
    double parsing_seconds = 0.0;
    double applicability_seconds = 0.0;
    double computation_seconds = 0.0;
    double writing_seconds = 0.0;
    double peak_resident_memory_mb = 0.0;
};

struct CalculationProvenance {
    RequestedCalculationProvenance requested;
    EffectiveCalculationProvenance effective;
    std::optional<ExecutionMetrics> execution_metrics;
};

struct ChargeResultRecord {
    MoleculeRecordIdentity identity;
    std::optional<charges::ChargeSet> charges;
    calculation::ExecutionStatus status = calculation::ExecutionStatus::success;
    std::vector<ResultDiagnostic> diagnostics;
};

struct ChargeResultDocument {
    std::string generator_name;
    std::string generator_version;
    calculation::ExecutionStatus status = calculation::ExecutionStatus::success;
    std::vector<ResultDiagnostic> diagnostics;
    std::vector<ChargeResultRecord> records;
    std::optional<CalculationProvenance> calculation_provenance;
};

[[nodiscard]] auto make_charge_result_document(
    std::span<const ImportedMoleculeRecord> records,
    const RequestedCalculationProvenance& requested, const calculation::ExecutionResult& result,
    std::string_view generator_name, std::string_view generator_version,
    std::optional<ExecutionMetrics> execution_metrics = std::nullopt) -> ChargeResultDocument;

} // namespace chargefw::adapters
