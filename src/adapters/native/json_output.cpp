#include <chargefw/adapters/native/json_output.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <string>

namespace chargefw::adapters::native::json_output {
namespace {

using Json = nlohmann::json;

constexpr auto charge_decimal_places = 4;
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

[[nodiscard]] auto mapping_json(const MoleculeRecordMapping& mapping, const bool atoms) -> Json {
    const auto& values = atoms ? mapping.atom_indices : mapping.conformer_indices;
    if (values.empty()) {
        return Json{{"kind", "identity"}};
    }

    const auto identity =
        std::ranges::all_of(values, [index = std::size_t{0}](const auto& value) mutable {
            return value.has_value() && *value == index++;
        });
    if (identity) {
        return Json{{"kind", "identity"}};
    }

    Json encoded_values = Json::array();
    for (const auto& value : values) {
        encoded_values.push_back(value.has_value() ? Json{*value} : Json{nullptr});
    }
    return Json{{"kind", "source_to_native"}, {"values", std::move(encoded_values)}};
}

[[nodiscard]] auto record_json(const ChargeResultRecord& record, const std::size_t molecule_index)
    -> Json {
    Json input{{"source", record.identity.source},
               {"record_index", record.identity.record_index},
               {"atom_mapping", mapping_json(record.mapping, true)},
               {"conformer_mapping", mapping_json(record.mapping, false)}};
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

    Json calculation{{"method", Json{{"id", record.charges->method_id()}}}};
    if (const auto parameter_set_id = record.charges->parameter_set_id();
        parameter_set_id.has_value()) {
        calculation["parameter_set"] = Json{{"id", *parameter_set_id}};
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
    result["calculation"] = std::move(calculation);
    result["assignments"] = std::move(assignments);
    result["diagnostics"] = Json::array();
    return result;
}

} // namespace

JsonWriter::JsonWriter(std::ostream& output) : output_{std::addressof(output)} {}

auto JsonWriter::write(const ChargeResultDocument& document) const -> void {
    Json records = Json::array();
    for (std::size_t index = 0; index < document.records.size(); ++index) {
        records.push_back(record_json(document.records[index], index));
    }

    *output_ << Json{{"schema_version", "1.0"},
                     {"generator",
                      {{"name", document.generator_name}, {"version", document.generator_version}}},
                     {"results", std::move(records)}}
                    .dump(2)
             << '\n';
}

} // namespace chargefw::adapters::native::json_output
