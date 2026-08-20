#include <chargefw/adapters/native/json_output.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstddef>
#include <ostream>
#include <print>
#include <stdexcept>
#include <string>

namespace chargefw::adapters::native::json_output {
namespace {

using Json = nlohmann::json;

constexpr auto charge_scale = 10000.0;

[[nodiscard]] auto rounded_charge(const double value) -> double {
    return std::round(value * charge_scale) / charge_scale;
}

[[nodiscard]] auto charges_json(const charges::AtomicCharges& charges) -> Json {
    Json values = Json::array();
    for (const auto value : charges.values()) {
        values.push_back(rounded_charge(value));
    }
    return values;
}

[[nodiscard]] auto record_json(const ChargeResultRecord& record, const std::size_t molecule_index)
    -> Json {
    Json input{{"source", record.identity.source}, {"record_index", record.identity.record_index}};
    if (!record.identity.record_id.empty()) {
        input["record_id"] = record.identity.record_id;
    }

    Json result{{"input", std::move(input)}};
    if (!record.charges.has_value()) {
        result["status"] = "error";
        result["diagnostics"] =
            Json::array({{{"severity", "error"},
                          {"code", "no_applicable_method"},
                          {"message", "No applicable method and parameter set were found."}}});
        return result;
    }

    Json assignments = Json::array();
    for (const auto& assignment : record.charges->assignments()) {
        if (assignment.target.molecule_index != molecule_index) {
            continue;
        }

        Json encoded_assignment{{"atom_order", "source"},
                                {"charges", charges_json(assignment.charges)},
                                {"total_charge", assignment.charges.total()}};
        if (assignment.target.conformer_index.has_value()) {
            encoded_assignment["conformer_index"] = *assignment.target.conformer_index;
        }
        assignments.push_back(std::move(encoded_assignment));
    }
    if (assignments.empty()) {
        throw std::invalid_argument{"Charge result record has no assignment for molecule index " +
                                    std::to_string(molecule_index)};
    }

    result["status"] = "success";
    result["assignments"] = std::move(assignments);
    result["diagnostics"] = Json::array();
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

    Json requested{
        {"method", optional_id(provenance.requested.method_id)},
        {"parameter_set", optional_id(provenance.requested.parameter_set_id)},
        {"classification", {{"permissive_types", provenance.requested.permissive_types}}},
        {"resource_policy",
         provenance.requested.full_atom_threshold.has_value()
             ? Json{{"full_atom_threshold", *provenance.requested.full_atom_threshold},
                    {"max_threads", provenance.requested.max_threads}}
             : Json{{"full_atom_threshold", "unlimited"},
                    {"max_threads", provenance.requested.max_threads}}},
        {"execution",
         {{"kind", provenance.requested.execution_kind},
          {"radius_angstrom", optional_value(provenance.requested.execution_radius)},
          {"charge_correction",
           optional_value(provenance.requested.execution_charge_correction)}}}};
    requested["input"] = {{"conformers", provenance.requested.conformer_selection}};
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
    return {{"requested", std::move(requested)}, {"effective", std::move(effective)}};
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
        {"results", std::move(records)}};
    if (document.calculation_provenance.has_value()) {
        result["calculation_provenance"] = provenance_json(*document.calculation_provenance);
    }
    std::print(*output_, "{}\n", result.dump(2));
}

} // namespace chargefw::adapters::native::json_output
