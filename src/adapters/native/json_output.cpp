#include <chargefw/adapters/native/json_output.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstddef>
#include <ostream>
#include <print>
#include <span>
#include <stdexcept>
#include <string>

namespace chargefw::adapters::native::json_output {
namespace {

using Json = nlohmann::json;

constexpr auto charge_scale = 10000.0;
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
        values.push_back(rounded(value, charge_scale));
    }
    return values;
}

[[nodiscard]] auto record_json(const ChargeResultRecord& record, const std::size_t molecule_index)
    -> Json {
    Json input{{"source", record.identity.source}, {"record_index", record.identity.record_index}};
    if (!record.identity.record_id.empty()) {
        input["record_id"] = record.identity.record_id;
    }

    Json result{{"input", std::move(input)}, {"status", calculation::to_string(record.status)}};
    if (!record.charges.has_value()) {
        result["diagnostics"] = diagnostics_json(record.diagnostics);
        return result;
    }

    Json assignments = Json::array();
    for (const auto& assignment : record.charges->assignments()) {
        if (assignment.target.molecule_index != molecule_index) {
            continue;
        }

        auto encoded_charges = charges_json(assignment.charges);
        auto serialized_total = 0.0;
        for (const auto& value : encoded_charges) {
            serialized_total += value.get<double>();
        }
        Json encoded_assignment{{"atom_order", "source"},
                                {"charges", std::move(encoded_charges)},
                                {"total_charge", rounded(serialized_total, charge_scale)}};
        if (assignment.target.conformer_index.has_value()) {
            encoded_assignment["conformer_index"] = *assignment.target.conformer_index;
        }
        assignments.push_back(std::move(encoded_assignment));
    }
    if (assignments.empty()) {
        throw std::invalid_argument{"Charge result record has no assignment for molecule index " +
                                    std::to_string(molecule_index)};
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
          {"radius_angstrom", optional_value(provenance.requested.execution_radius)},
          {"charge_correction",
           optional_value(provenance.requested.execution_charge_correction)}}}};
    if (provenance.requested.conformer_selection.has_value()) {
        requested["input"] = {{"conformers", *provenance.requested.conformer_selection}};
    }
    requested["method_options"] = method_options_json(provenance.requested.method_options);
    if (provenance.requested.structural_input_policy.has_value()) {
        requested["structural_input"] = {
            {"selection", provenance.requested.structural_input_policy->selection},
            {"bonds", provenance.requested.structural_input_policy->bonds}};
    }

    Json effective{{"method", optional_id(provenance.effective.method_id)},
                   {"parameter_set", optional_id(provenance.effective.parameter_set_id)},
                   {"warnings", provenance.effective.warnings}};
    if (provenance.effective.execution_mode.has_value()) {
        effective["execution"] = {
            {"mode", *provenance.effective.execution_mode},
            {"radius_angstrom", optional_value(provenance.effective.execution_radius)},
            {"charge_correction",
             optional_value(provenance.effective.execution_charge_correction)}};
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
    for (std::size_t index = 0; index < document.records.size(); ++index) {
        records.push_back(record_json(document.records[index], index));
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

} // namespace chargefw::adapters::native::json_output
