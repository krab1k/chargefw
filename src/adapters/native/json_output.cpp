#include <chargefw/adapters/native/json_output.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstddef>
#include <ostream>
#include <print>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::adapters::native::json_output {
namespace {

using Json = nlohmann::json;

constexpr auto metric_scale = 1000.0;

[[nodiscard]] auto diagnostics_json(const std::span<const ResultDiagnostic> diagnostics) -> Json {
    const auto severity_name = [](const DiagnosticSeverity severity) -> const char* {
        switch (severity) {
        case DiagnosticSeverity::info:
            return "info";
        case DiagnosticSeverity::warning:
            return "warning";
        case DiagnosticSeverity::error:
            return "error";
        }
        throw std::invalid_argument{"unknown diagnostic severity"};
    };
    Json result = Json::array();
    for (const auto& diagnostic : diagnostics) {
        auto encoded = Json{{"severity", severity_name(diagnostic.severity)},
                            {"code", diagnostic.code},
                            {"message", diagnostic.message}};
        const auto add_index = [&encoded](const std::string_view name,
                                          const std::optional<std::size_t> index) {
            if (index.has_value()) {
                encoded[name] = *index;
            }
        };
        add_index("molecule_index", diagnostic.molecule_index);
        add_index("atom_index", diagnostic.atom_index);
        add_index("bond_index", diagnostic.bond_index);
        add_index("conformer_index", diagnostic.conformer_index);
        add_index("line", diagnostic.line);
        result.push_back(std::move(encoded));
    }
    return result;
}

[[nodiscard]] auto rounded(const double value, const double scale) -> double {
    return std::round(value * scale) / scale;
}

[[nodiscard]] auto charges_json(const charges::AtomicCharges& charges) -> Json {
    Json values = Json::array();
    for (const auto value : charges.values()) {
        values.push_back(value);
    }
    return values;
}

[[nodiscard]] auto source_format_name(const MolecularSourceFormat format) -> std::string_view {
    switch (format) {
    case MolecularSourceFormat::molecule_json:
        return "molecule-json";
    case MolecularSourceFormat::mol:
        return "mol";
    case MolecularSourceFormat::sdf:
        return "sdf";
    case MolecularSourceFormat::mol2:
        return "mol2";
    case MolecularSourceFormat::pdb:
        return "pdb";
    case MolecularSourceFormat::mmcif:
        return "mmcif";
    }
    throw std::invalid_argument{"unknown molecular source format"};
}

[[nodiscard]] auto connectivity_name(const SourceConnectivity connectivity) -> std::string_view {
    switch (connectivity) {
    case SourceConnectivity::absent:
        return "absent";
    case SourceConnectivity::explicitly_empty:
        return "explicitly-empty";
    case SourceConnectivity::present:
        return "present";
    }
    throw std::invalid_argument{"unknown source connectivity state"};
}

[[nodiscard]] auto hierarchy_json(const SourceHierarchyLabels& labels) -> Json {
    auto result = Json::object();
    const auto add = [&result](const std::string_view name,
                               const std::optional<std::string>& value) {
        if (value.has_value()) {
            result[name] = *value;
        }
    };
    add("atom", labels.atom);
    add("residue", labels.residue);
    add("chain", labels.chain);
    add("sequence", labels.sequence);
    return result;
}

[[nodiscard]] auto source_reference_json(const SourceAtomReference& reference) -> Json {
    auto result = Json{{"source_position", reference.position}};
    if (reference.id.has_value()) {
        result["source_id"] = *reference.id;
    }
    if (reference.structural_labels.has_value()) {
        const auto& labels = *reference.structural_labels;
        auto structural = Json{{"author", hierarchy_json(labels.author)},
                               {"label", hierarchy_json(labels.label)}};
        const auto add = [&structural](const std::string_view name,
                                       const std::optional<std::string>& value) {
            if (value.has_value()) {
                structural[name] = *value;
            }
        };
        add("entity", labels.entity);
        add("insertion_code", labels.insertion_code);
        add("alternate_location", labels.alternate_location);
        add("segment", labels.segment);
        result["structural_labels"] = std::move(structural);
    }
    return result;
}

[[nodiscard]] auto import_metadata_json(const MoleculeImportMetadata& metadata) -> Json {
    auto atoms = Json::array();
    for (const auto& reference : metadata.atoms) {
        atoms.push_back(source_reference_json(reference));
    }
    auto conformers = Json::array();
    for (const auto& conformer : metadata.conformers) {
        auto sites = Json::array();
        for (const auto& reference : conformer.sites) {
            sites.push_back(source_reference_json(reference));
        }
        auto encoded = Json{{"source_position", conformer.position}, {"sites", std::move(sites)}};
        if (conformer.id.has_value()) {
            encoded["source_id"] = *conformer.id;
        }
        conformers.push_back(std::move(encoded));
    }

    auto policy = Json::object();
    const auto add_policy = [&policy](const std::string_view name,
                                      const std::optional<std::string>& value) {
        if (value.has_value()) {
            policy[name] = *value;
        }
    };
    add_policy("record_selection", metadata.record_selection);
    add_policy("alternate_location_selection", metadata.alternate_location_selection);
    add_policy("conformer_selection", metadata.conformer_selection);
    add_policy("bond_strategy", metadata.bond_strategy);

    return Json{{"format", source_format_name(metadata.format)},
                {"policy", std::move(policy)},
                {"source_connectivity", connectivity_name(metadata.source_connectivity)},
                {"atom_mapping", std::move(atoms)},
                {"conformer_mapping", std::move(conformers)}};
}

[[nodiscard]] auto record_json(const ChargeResultRecord& record) -> Json {
    Json input{{"source", record.input.identity.source},
               {"record_index", record.input.identity.record_index}};
    if (!record.input.identity.record_id.empty()) {
        input["record_id"] = record.input.identity.record_id;
    }
    if (record.input.import_metadata.has_value()) {
        input["import"] = import_metadata_json(*record.input.import_metadata);
    }

    Json result{{"input", std::move(input)}, {"status", calculation::to_string(record.status)}};
    if (record.assignments.empty()) {
        result["diagnostics"] = diagnostics_json(record.diagnostics);
        return result;
    }

    Json assignments = Json::array();
    for (const auto& assignment : record.assignments) {
        auto encoded_charges = charges_json(assignment.charges);
        auto total_charge = 0.0;
        for (const auto value : assignment.charges.values()) {
            total_charge += value;
        }
        Json target{{"molecule_index", assignment.target.molecule_index}};
        if (assignment.target.conformer_index.has_value()) {
            target["conformer_index"] = *assignment.target.conformer_index;
        }
        Json encoded_assignment{
            {"scope", assignment.target.conformer_index.has_value() ? "conformer" : "molecule"},
            {"target", std::move(target)},
            {"charge_unit", "e"},
            {"charges", std::move(encoded_charges)},
            {"total_charge", total_charge}};
        assignments.push_back(std::move(encoded_assignment));
    }
    result["assignments"] = std::move(assignments);
    result["diagnostics"] = diagnostics_json(record.diagnostics);
    return result;
}

[[nodiscard]] auto provenance_json(const CalculationProvenance& provenance) -> Json {
    const auto optional_id = [](const std::optional<std::string>& id) -> Json {
        return id.has_value() ? Json{{"id", *id}} : Json(nullptr);
    };
    const auto optional_value = [](const auto& value) -> Json {
        return value.has_value() ? Json(*value) : Json(nullptr);
    };
    const auto option_value = [](const methods::MethodOptionValue& value) -> Json {
        return std::visit([](const auto& item) -> Json { return item; }, value);
    };
    const auto method_options_json =
        [&option_value](const std::map<std::string, methods::MethodOptions>& options) -> Json {
        Json result = Json::object();
        for (const auto& [method_id, values] : options) {
            Json method = Json::object();
            for (const auto& [id, value] : values.values()) {
                method[id] = option_value(value);
            }
            result[method_id] = std::move(method);
        }
        return result;
    };
    const auto threshold_value = [](const std::optional<std::size_t>& threshold) -> Json {
        return threshold.has_value() ? Json(*threshold) : Json("unlimited");
    };

    Json requested{
        {"method", optional_id(provenance.requested.method_id)},
        {"parameter_set", optional_id(provenance.requested.parameter_set_id)},
        {"classification", {{"permissive_types", provenance.requested.permissive_types}}},
        {"resource_policy",
         {{"cutoff_atom_threshold", threshold_value(provenance.requested.cutoff_atom_threshold)},
          {"cover_atom_threshold", threshold_value(provenance.requested.cover_atom_threshold)},
          {"max_threads", provenance.requested.max_threads}}},
        {"execution",
         {{"kind", provenance.requested.execution_kind},
          {"radius_angstrom", optional_value(provenance.requested.execution_radius)}}}};
    requested["method_options"] = method_options_json(provenance.requested.method_options);

    Json effective{{"method", optional_id(provenance.effective.method_id)},
                   {"parameter_set", optional_id(provenance.effective.parameter_set_id)},
                   {"warnings", provenance.effective.warnings}};
    if (provenance.effective.execution_mode.has_value()) {
        effective["execution"] = {
            {"mode", *provenance.effective.execution_mode},
            {"radius_angstrom", optional_value(provenance.effective.execution_radius)}};
    }
    effective["method_options"] = method_options_json(provenance.effective.method_options);
    Json result{{"requested", std::move(requested)}, {"effective", std::move(effective)}};
    if (provenance.execution_metrics.has_value()) {
        const auto& metrics = *provenance.execution_metrics;
        result["execution_metrics"] = {
            {"started_at", metrics.started_at},
            {"ended_at", metrics.ended_at},
            {"runtime_seconds", rounded(metrics.runtime_seconds, metric_scale)},
            {"phases",
             {{"parsing_seconds", rounded(metrics.parsing_seconds, metric_scale)},
              {"applicability_seconds", rounded(metrics.applicability_seconds, metric_scale)},
              {"computation_seconds", rounded(metrics.computation_seconds, metric_scale)},
              {"writing_seconds", rounded(metrics.writing_seconds, metric_scale)}}},
            {"peak_resident_memory_mb", rounded(metrics.peak_resident_memory_mb, metric_scale)}};
    }
    return result;
}

} // namespace

JsonWriter::JsonWriter(std::ostream& output) : output_{std::addressof(output)} {}

auto JsonWriter::write(const ChargeResultDocument& document) const -> void {
    Json records = Json::array();
    for (const auto& record : document.records) {
        records.push_back(record_json(record));
    }

    Json result{
        {"schema_version", "1.0"},
        {"generator", {{"name", document.generator_name}, {"version", document.generator_version}}},
        {"status", calculation::to_string(document.status)},
        {"diagnostics", diagnostics_json(document.diagnostics)},
        {"results", std::move(records)}};
    if (document.calculation_provenance.has_value()) {
        result["calculation_provenance"] = provenance_json(*document.calculation_provenance);
    }
    std::print(*output_, "{}\n", result.dump(2));
}

auto JsonWriter::write(const ChargeCalculationResult& result, const std::string_view generator_name,
                       const std::string_view generator_version,
                       std::optional<ExecutionMetrics> execution_metrics) const -> void {
    write(make_charge_result_document(result, generator_name, generator_version,
                                      std::move(execution_metrics)));
}

} // namespace chargefw::adapters::native::json_output
