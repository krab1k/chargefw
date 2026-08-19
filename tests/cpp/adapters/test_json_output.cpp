#include <cassert>
#include <chargefw/adapters/charge_result_document.h>
#include <chargefw/adapters/native/json_output.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>

#include <nlohmann/json.hpp>

#include <sstream>
#include <string>
#include <vector>

namespace adapters = chargefw::adapters;
namespace charges = chargefw::charges;
namespace json_output = chargefw::adapters::native::json_output;

auto main() -> int {
    const auto document = adapters::ChargeResultDocument{
        .generator_name = "ChargeFW",
        .generator_version = "test",
        .records =
            {{.identity = {.source = "water.sdf", .record_index = 2, .record_id = "water"},
              .mapping = {.atom_indices = {0, 1, 2}, .conformer_indices = {0}},
              .charges = charges::ChargeSet{"formal",
                                            {{.target = {.molecule_index = 0, .conformer_index = 0},
                                              .charges = charges::AtomicCharges{{-0.87654, 0.43827,
                                                                                 0.43827}}}}}},
             {.identity = {.source = "water.sdf", .record_index = 3, .record_id = "unavailable"},
              .mapping = {.atom_indices = {}, .conformer_indices = {}},
              .charges = std::nullopt}},
        .calculation_provenance = adapters::CalculationProvenance{
            .effective_execution_mode = "cutoff",
            .radius = 8.0,
            .charge_correction = "uniform",
            .permissive_types = true,
            .full_atom_threshold = std::nullopt,
            .warnings = {"full execution exceeds the shared threshold"}}};

    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(document);
    const auto result = nlohmann::json::parse(output.str());

    assert(result.at("schema_version") == "1.0");
    assert(result.at("generator").at("name") == "ChargeFW");
    assert(result.at("results").size() == 2);
    const auto& provenance = result.at("calculation_provenance");
    assert(provenance.at("execution").at("mode") == "cutoff");
    assert(provenance.at("execution").at("radius_angstrom") == 8.0);
    assert(provenance.at("execution").at("charge_correction") == "uniform");
    assert(provenance.at("classification").at("permissive_types") == true);
    assert(provenance.at("resource_policy").at("full_atom_threshold") == "unlimited");
    assert(provenance.at("warnings").at(0) == "full execution exceeds the shared threshold");

    const auto& calculated = result.at("results").at(0);
    assert(calculated.at("status") == "success");
    assert(calculated.at("input").at("atom_mapping").at("kind") == "identity");
    assert(calculated.at("calculation").at("method").at("id") == "formal");
    assert(calculated.at("assignments").at(0).at("conformer_index") == 0);
    assert(calculated.at("assignments").at(0).at("charges").size() == 3);
    assert(calculated.at("assignments").at(0).at("charges").at(0) == -0.8765);
    assert(calculated.at("assignments").at(0).at("charges").at(1) == 0.4383);

    const auto& unavailable = result.at("results").at(1);
    assert(unavailable.at("status") == "error");
    assert(unavailable.at("input").at("atom_mapping").at("kind") == "identity");
    assert(unavailable.at("diagnostics").at(0).at("code") == "no_applicable_method");

    return 0;
}
