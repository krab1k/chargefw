#include <chargefw/adapters/charge_result_document.h>
#include <chargefw/adapters/native/json_output.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <snitch/snitch.hpp>

#include <nlohmann/json.hpp>

#include <sstream>
#include <string>
#include <vector>

namespace adapters = chargefw::adapters;
namespace calculation = chargefw::calculation;
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
              .charges = std::nullopt,
              .status = calculation::ExecutionStatus::no_executable_plan,
              .diagnostics = {{.severity = adapters::DiagnosticSeverity::error,
                               .code = "no_executable_plan",
                               .message = "No executable plan.",
                               .molecule_index = 1,
                               .atom_index = 4}}}},
        .calculation_provenance = adapters::CalculationProvenance{
            .requested = {.method_id = std::nullopt,
                          .parameter_set_id = std::nullopt,
                          .permissive_types = true,
                          .cutoff_atom_threshold = std::nullopt,
                          .cover_atom_threshold = std::nullopt,
                          .max_threads = 0,
                          .execution_kind = "auto",
                          .execution_radius = std::nullopt,
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
    CHECK(result.at("status") == "success");
    CHECK(result.at("diagnostics").empty());
    CHECK(result.at("generator").at("name") == "ChargeFW");
    REQUIRE(result.at("results").size() == 2);
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
    CHECK(effective.at("warnings").at(0) == "full execution exceeds the shared threshold");
    CHECK(effective.at("method_options").at("formal").empty());

    const auto& calculated = result.at("results").at(0);
    CHECK(calculated.at("status") == "success");
    CHECK_FALSE(calculated.at("input").contains("atom_mapping"));
    const auto& assignment = calculated.at("assignments").at(0);
    CHECK(assignment.at("scope") == "conformer");
    CHECK(assignment.at("target").at("molecule_index") == 0);
    CHECK(assignment.at("target").at("conformer_index") == 0);
    CHECK(assignment.at("charge_unit") == "e");
    REQUIRE(assignment.at("charges").size() == 3);
    CHECK(assignment.at("charges").at(0) == -0.87654);
    CHECK(assignment.at("charges").at(1) == 0.43827);
    CHECK(assignment.at("total_charge") == 0.0);

    const auto& unavailable = result.at("results").at(1);
    CHECK(unavailable.at("status") == "no_executable_plan");
    CHECK_FALSE(unavailable.at("input").contains("atom_mapping"));
    CHECK(unavailable.at("diagnostics").at(0).at("code") == "no_executable_plan");
    CHECK(unavailable.at("diagnostics").at(0).at("molecule_index") == 1);
    CHECK(unavailable.at("diagnostics").at(0).at("atom_index") == 4);
}

TEST_CASE("JSON output serializes a cancelled result without assignments", "[adapters][json]") {
    const auto diagnostic =
        adapters::ResultDiagnostic{.severity = adapters::DiagnosticSeverity::info,
                                   .code = "calculation_cancelled",
                                   .message = "Calculation was cancelled."};
    const auto document = adapters::ChargeResultDocument{
        .generator_name = "ChargeFW",
        .generator_version = "test",
        .status = calculation::ExecutionStatus::cancelled,
        .diagnostics = {diagnostic},
        .records = {{.identity = {.source = "water.sdf", .record_index = 0},
                     .charges = std::nullopt,
                     .status = calculation::ExecutionStatus::cancelled,
                     .diagnostics = {diagnostic}}}};

    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(document);
    const auto result = nlohmann::json::parse(output.str());

    CHECK(result.at("status") == "cancelled");
    CHECK(result.at("diagnostics").at(0).at("code") == "calculation_cancelled");
    const auto& record = result.at("results").at(0);
    CHECK(record.at("status") == "cancelled");
    CHECK_FALSE(record.contains("assignments"));
}

TEST_CASE("result assembly validates assignment dimensions targets and scope", "[adapters][json]") {
    const auto records = std::vector{adapters::ImportedMoleculeRecord{
        .molecule =
            chargefw::core::Molecule{
                std::vector{chargefw::core::Atom{1}, chargefw::core::Atom{1}},
                {},
                std::vector{chargefw::core::Conformer{{chargefw::core::Position{0.0, 0.0, 0.0},
                                                       chargefw::core::Position{1.0, 0.0, 0.0}}}},
                "hydrogen"},
        .identity = {.source = "hydrogen.json", .record_id = "hydrogen"}}};
    const auto make_document = [&records](calculation::ExecutionResult result) {
        if (result.status == calculation::ExecutionStatus::success &&
            !result.effective.has_value()) {
            result.effective = calculation::EffectiveCalculation{
                .method_id = "formal", .execution_policy = calculation::ExecutionPolicy{}};
        }
        return adapters::make_charge_result_document(records, {}, result, "ChargeFW", "test");
    };

    CHECK_THROWS_AS(make_document(calculation::ExecutionResult{}), std::invalid_argument);
    CHECK_THROWS_AS(
        make_document(calculation::ExecutionResult{
            .charges = charges::ChargeSet{"formal",
                                          {{.target = {.molecule_index = 0},
                                            .charges = charges::AtomicCharges{{0.0}}}}}}),
        std::invalid_argument);
    CHECK_THROWS_AS(
        make_document(calculation::ExecutionResult{
            .charges = charges::ChargeSet{"formal",
                                          {{.target = {.molecule_index = 1},
                                            .charges = charges::AtomicCharges{{0.0, 0.0}}}}}}),
        std::invalid_argument);
    CHECK_THROWS_AS(
        make_document(calculation::ExecutionResult{
            .charges = charges::ChargeSet{"formal",
                                          {{.target = {.molecule_index = 0, .conformer_index = 1},
                                            .charges = charges::AtomicCharges{{0.0, 0.0}}}}}}),
        std::invalid_argument);
    CHECK_THROWS_AS(
        make_document(calculation::ExecutionResult{
            .charges = charges::ChargeSet{"formal",
                                          {{.target = {.molecule_index = 0},
                                            .charges = charges::AtomicCharges{{0.0, 0.0}}},
                                           {.target = {.molecule_index = 0, .conformer_index = 0},
                                            .charges = charges::AtomicCharges{{0.0, 0.0}}}}}}),
        std::invalid_argument);
    CHECK_THROWS_AS(
        make_document(calculation::ExecutionResult{
            .status = calculation::ExecutionStatus::numerical_failure,
            .charges = charges::ChargeSet{"formal",
                                          {{.target = {.molecule_index = 0},
                                            .charges = charges::AtomicCharges{{0.0, 0.0}}}}}}),
        std::invalid_argument);
}
