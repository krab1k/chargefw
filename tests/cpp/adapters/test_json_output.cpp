#include <chargefw/adapters/charge_result_document.h>
#include <chargefw/adapters/native/json_output.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>
#include <snitch/snitch.hpp>

#include <nlohmann/json.hpp>

#include <sstream>
#include <string>
#include <vector>

namespace adapters = chargefw::adapters;
namespace charges = chargefw::charges;
namespace json_output = chargefw::adapters::native::json_output;

TEST_CASE("JSON output serializes ordered records and calculation provenance", "[adapters][json]") {
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
                          .cutoff_atom_threshold = std::nullopt,
                          .cover_atom_threshold = std::nullopt,
                          .max_threads = 0,
                          .execution_kind = "auto",
                          .execution_radius = std::nullopt,
                          .execution_charge_correction = std::nullopt,
                          .structural_input_policy =
                              adapters::StructuralInputPolicyProvenance{.selection = "polymers",
                                                                        .bonds = "hybrid"},
                          .conformer_selection = "all",
                          .method_options = {{"qeq",
                                              chargefw::methods::MethodOptions{
                                                  {{"overlap_term", "Ohno"}}}}}},
            .effective = {.method_id = "formal",
                          .parameter_set_id = "test-formal",
                          .execution_mode = "cutoff",
                          .execution_radius = 8.0,
                          .execution_charge_correction = "uniform",
                          .warnings = {"full execution exceeds the shared threshold"},
                          .method_options = {{"formal", chargefw::methods::MethodOptions{}}}},
            .execution_metrics =
                adapters::ExecutionMetrics{.started_at = "2026-08-20T10:00:00.000Z",
                                           .ended_at = "2026-08-20T10:00:01.250Z",
                                           .runtime_seconds = 1.23456,
                                           .parsing_seconds = 0.1004,
                                           .applicability_seconds = 0.20,
                                           .computation_seconds = 0.80,
                                           .writing_seconds = 0.15,
                                           .peak_resident_memory_mb = 123.4567}}};

    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(document);
    const auto result = nlohmann::json::parse(output.str());

    CHECK(result.at("schema_version") == "1.0");
    CHECK(result.at("generator").at("name") == "ChargeFW");
    CHECK(result.at("results").size() == 2);
    const auto& provenance = result.at("calculation_provenance");
    const auto& requested = provenance.at("requested");
    CHECK(requested.at("method").is_null());
    CHECK(requested.at("parameter_set").is_null());
    CHECK(requested.at("execution").at("kind") == "auto");
    CHECK(requested.at("execution").at("radius_angstrom").is_null());
    CHECK(requested.at("classification").at("permissive_types") == true);
    CHECK(requested.at("resource_policy").at("cutoff_atom_threshold") == "unlimited");
    CHECK(requested.at("resource_policy").at("cover_atom_threshold") == "unlimited");
    CHECK(requested.at("resource_policy").at("max_threads") == 0);
    CHECK(requested.at("structural_input").at("selection") == "polymers");
    CHECK(requested.at("structural_input").at("bonds") == "hybrid");
    CHECK(requested.at("input").at("conformers") == "all");
    CHECK(requested.at("method_options").at("qeq").at("overlap_term") == "Ohno");
    const auto& effective = provenance.at("effective");
    const auto& metrics = provenance.at("execution_metrics");
    CHECK(metrics.at("started_at") == "2026-08-20T10:00:00.000Z");
    CHECK(metrics.at("ended_at") == "2026-08-20T10:00:01.250Z");
    CHECK(metrics.at("runtime_seconds") == 1.235);
    CHECK(metrics.at("phases").at("parsing_seconds") == 0.1);
    CHECK(metrics.at("phases").at("applicability_seconds") == 0.20);
    CHECK(metrics.at("phases").at("computation_seconds") == 0.80);
    CHECK(metrics.at("phases").at("writing_seconds") == 0.15);
    CHECK(metrics.at("peak_resident_memory_mb") == 123.457);
    CHECK(effective.at("execution").at("mode") == "cutoff");
    CHECK(effective.at("method").at("id") == "formal");
    CHECK(effective.at("parameter_set").at("id") == "test-formal");
    CHECK(effective.at("execution").at("radius_angstrom") == 8.0);
    CHECK(effective.at("execution").at("charge_correction") == "uniform");
    CHECK(effective.at("warnings").at(0) == "full execution exceeds the shared threshold");
    CHECK(effective.at("method_options").at("formal").empty());

    const auto& calculated = result.at("results").at(0);
    CHECK(calculated.at("status") == "success");
    CHECK_FALSE(calculated.at("input").contains("atom_mapping"));
    CHECK(calculated.at("assignments").at(0).at("conformer_index") == 0);
    CHECK(calculated.at("assignments").at(0).at("charges").size() == 3);
    CHECK(calculated.at("assignments").at(0).at("charges").at(0) == -0.8765);
    CHECK(calculated.at("assignments").at(0).at("charges").at(1) == 0.4383);

    const auto& unavailable = result.at("results").at(1);
    CHECK(unavailable.at("status") == "error");
    CHECK_FALSE(unavailable.at("input").contains("atom_mapping"));
    CHECK(unavailable.at("diagnostics").at(0).at("code") == "no_applicable_method");
}
