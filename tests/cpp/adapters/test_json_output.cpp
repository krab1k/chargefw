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
              .charges = charges::ChargeSet{"formal",
                                            {{.target = {.molecule_index = 0, .conformer_index = 0},
                                              .charges = charges::AtomicCharges{{-0.87654, 0.43827,
                                                                                 0.43827}}}}}},
             {.identity = {.source = "water.sdf", .record_index = 3, .record_id = "unavailable"},
              .charges = std::nullopt}},
        .calculation_provenance = adapters::CalculationProvenance{
            .requested = {.method_id = std::nullopt,
                          .parameter_set_id = std::nullopt,
                          .permissive_types = true,
                          .full_atom_threshold = std::nullopt,
                          .execution_kind = "auto",
                          .execution_radius = std::nullopt,
                          .execution_charge_correction = std::nullopt,
                          .structural_input_policy =
                              adapters::StructuralInputPolicyProvenance{.selection = "polymers",
                                                                        .bonds = "hybrid"},
                          .conformer_selection = "all"},
            .effective = {.method_id = "formal",
                          .parameter_set_id = "test-formal",
                          .execution_mode = "cutoff",
                          .execution_radius = 8.0,
                          .execution_charge_correction = "uniform",
                          .warnings = {"full execution exceeds the shared threshold"}}}};

    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(document);
    const auto result = nlohmann::json::parse(output.str());

    assert(result.at("schema_version") == "1.0");
    assert(result.at("generator").at("name") == "ChargeFW");
    assert(result.at("results").size() == 2);
    const auto& provenance = result.at("calculation_provenance");
    const auto& requested = provenance.at("requested");
    assert(requested.at("method").is_null());
    assert(requested.at("parameter_set").is_null());
    assert(requested.at("execution").at("kind") == "auto");
    assert(requested.at("execution").at("radius_angstrom").is_null());
    assert(requested.at("classification").at("permissive_types") == true);
    assert(requested.at("resource_policy").at("full_atom_threshold") == "unlimited");
    assert(requested.at("structural_input").at("selection") == "polymers");
    assert(requested.at("structural_input").at("bonds") == "hybrid");
    assert(requested.at("input").at("conformers") == "all");
    const auto& effective = provenance.at("effective");
    assert(effective.at("execution").at("mode") == "cutoff");
    assert(effective.at("method").at("id") == "formal");
    assert(effective.at("parameter_set").at("id") == "test-formal");
    assert(effective.at("execution").at("radius_angstrom") == 8.0);
    assert(effective.at("execution").at("charge_correction") == "uniform");
    assert(effective.at("warnings").at(0) == "full execution exceeds the shared threshold");

    const auto& calculated = result.at("results").at(0);
    assert(calculated.at("status") == "success");
    assert(!calculated.at("input").contains("atom_mapping"));
    assert(calculated.at("assignments").at(0).at("conformer_index") == 0);
    assert(calculated.at("assignments").at(0).at("charges").size() == 3);
    assert(calculated.at("assignments").at(0).at("charges").at(0) == -0.8765);
    assert(calculated.at("assignments").at(0).at("charges").at(1) == 0.4383);

    const auto& unavailable = result.at("results").at(1);
    assert(unavailable.at("status") == "error");
    assert(!unavailable.at("input").contains("atom_mapping"));
    assert(unavailable.at("diagnostics").at(0).at("code") == "no_applicable_method");

    return 0;
}
