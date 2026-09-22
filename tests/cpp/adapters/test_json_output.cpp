#include <chargefw/adapters/charge_result.h>
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
#include <utility>
#include <vector>

namespace adapters = chargefw::adapters;
namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace json_output = chargefw::adapters::native::json_output;

TEST_CASE("JSON output serializes ordered records and calculation provenance", "[adapters][json]") {
    auto sites = std::vector<adapters::SourceAtomReference>{
        {.position = 0, .id = "10"}, {.position = 1, .id = "20"}, {.position = 2, .id = "30"}};
    const auto owned = adapters::make_charge_calculation_result(
        {{.molecule =
              chargefw::core::Molecule{
                  std::vector{chargefw::core::Atom{8}, chargefw::core::Atom{1},
                              chargefw::core::Atom{1}},
                  {},
                  {chargefw::core::Conformer{{chargefw::core::Position{0.0, 0.0, 0.0},
                                              chargefw::core::Position{1.0, 0.0, 0.0},
                                              chargefw::core::Position{0.0, 1.0, 0.0}}}}},
          .identity = {.source = "water.sdf", .record_index = 2, .record_id = 42},
          .import_metadata =
              adapters::MoleculeImportMetadata{.format = adapters::MolecularSourceFormat::sdf,
                                               .atoms = sites,
                                               .conformers = {{.position = 0,
                                                               .sites = std::move(sites)}},
                                               .conformer_selection = "all",
                                               .source_connectivity =
                                                   adapters::SourceConnectivity::present}}},
        {.method_id = std::nullopt,
         .parameter_set_id = std::nullopt,
         .permissive_types = true,
         .cutoff_atom_threshold = std::nullopt,
         .cover_atom_threshold = std::nullopt,
         .max_threads = 0,
         .execution_kind = "auto",
         .execution_radius = std::nullopt,
         .method_options = {{"qeq", chargefw::methods::MethodOptions{{{"overlap_term", "Ohno"}}}}}},
        {.charges =
             charges::ChargeSet{"formal",
                                {{.target = {.molecule_index = 0, .conformer_index = 0},
                                  .charges = charges::AtomicCharges{{-0.87654, 0.43827, 0.43827}}}},
                                "test-formal"},
         .effective = calculation::EffectiveCalculation{
             .method_id = "formal",
             .parameter_set_id = "test-formal",
             .method_options = {},
             .execution_policy =
                 calculation::ExecutionPolicy{calculation::ExecutionMode::cutoff, 8.0},
             .execution_issues = {
                 {.kind = chargefw::methods::ExecutionIssueKind::resource_threshold_exceeded,
                  .message = "full execution exceeds the shared threshold"}}}});
    const auto execution_metrics =
        adapters::ExecutionMetrics{.started_at = "2026-08-20T10:00:00.000Z",
                                   .ended_at = "2026-08-20T10:00:01.250Z",
                                   .runtime_seconds = 1.23456,
                                   .parsing_seconds = 0.1004,
                                   .applicability_seconds = 0.20,
                                   .computation_seconds = 0.80,
                                   .peak_resident_memory_mb = 123.4567};

    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(owned, "test", execution_metrics);
    const auto result = nlohmann::json::parse(output.str());

    CHECK(result.at("schema_version") == "1.0");
    CHECK(result.at("status") == "success");
    CHECK(result.at("diagnostics").empty());
    CHECK(result.at("generator").at("name") == "ChargeFW");
    REQUIRE(result.at("results").size() == 1);
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
    CHECK_FALSE(requested.contains("structural_input"));
    CHECK_FALSE(requested.contains("input"));
    CHECK(requested.at("method_options").at("qeq").at("overlap_term") == "Ohno");
    const auto& effective = provenance.at("effective");
    const auto& encoded_metrics = provenance.at("execution_metrics");
    CHECK(encoded_metrics.at("started_at") == "2026-08-20T10:00:00.000Z");
    CHECK(encoded_metrics.at("ended_at") == "2026-08-20T10:00:01.250Z");
    CHECK(encoded_metrics.at("runtime_seconds") == 1.235);
    CHECK(encoded_metrics.at("phases").at("parsing_seconds") == 0.1);
    CHECK(encoded_metrics.at("phases").at("applicability_seconds") == 0.20);
    CHECK(encoded_metrics.at("phases").at("computation_seconds") == 0.80);
    CHECK(encoded_metrics.at("peak_resident_memory_mb") == 123.457);
    CHECK(effective.at("execution").at("mode") == "cutoff");
    CHECK(effective.at("method").at("id") == "formal");
    CHECK(effective.at("parameter_set").at("id") == "test-formal");
    CHECK(effective.at("execution").at("radius_angstrom") == 8.0);
    CHECK(effective.at("warnings").at(0) == "full execution exceeds the shared threshold");
    CHECK(effective.at("method_options").at("formal").empty());

    const auto& calculated = result.at("results").at(0);
    CHECK(calculated.at("status") == "success");
    CHECK(calculated.at("input").at("record_id") == 42);
    const auto& imported = calculated.at("input").at("import");
    CHECK(imported.at("format") == "sdf");
    CHECK(imported.at("policy").at("conformer_selection") == "all");
    CHECK(imported.at("source_connectivity") == "present");
    CHECK(imported.at("atom_mapping").at(0).at("source_id") == "10");
    CHECK(imported.at("atom_mapping").at(2).at("source_position") == 2);
    const auto& assignment = calculated.at("assignments").at(0);
    CHECK(assignment.at("scope") == "conformer");
    CHECK(assignment.at("target").at("molecule_index") == 0);
    CHECK(assignment.at("target").at("conformer_index") == 0);
    CHECK(assignment.at("charge_unit") == "e");
    REQUIRE(assignment.at("charges").size() == 3);
    CHECK(assignment.at("charges").at(0) == -0.87654);
    CHECK(assignment.at("charges").at(1) == 0.43827);
    CHECK(assignment.at("total_charge") == 0.0);
}

TEST_CASE("JSON output serializes a cancelled result without assignments", "[adapters][json]") {
    const auto owned = adapters::make_charge_calculation_result(
        {{.molecule = chargefw::core::Molecule{std::vector{chargefw::core::Atom{8}}},
          .identity = {.source = "water.sdf", .record_index = 0}}},
        {}, {.status = calculation::ExecutionStatus::cancelled});

    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(owned, "test");
    const auto result = nlohmann::json::parse(output.str());

    CHECK(result.at("status") == "cancelled");
    CHECK(result.at("diagnostics").at(0).at("code") == "calculation_cancelled");
    const auto& record = result.at("results").at(0);
    CHECK(record.at("status") == "cancelled");
    CHECK_FALSE(record.contains("assignments"));
}

TEST_CASE("JSON output serializes caller atom IDs for a manual record", "[adapters][json]") {
    const auto owned = adapters::make_charge_calculation_result(
        {{.molecule = chargefw::core::Molecule{std::vector{chargefw::core::Atom{1},
                                                           chargefw::core::Atom{1}},
                                               {},
                                               {},
                                               "hydrogen"},
          .identity = {.source = "manual", .record_index = 1},
          .caller_atom_ids = std::vector<adapters::PortableId>{"H", std::int64_t{9}}}},
        {},
        {.charges = charges::ChargeSet{"formal",
                                       {{.target = {.molecule_index = 0},
                                         .charges = charges::AtomicCharges{{0.0, 0.0}}}},
                                       std::nullopt},
         .effective = calculation::EffectiveCalculation{
             .method_id = "formal", .execution_policy = calculation::ExecutionPolicy{}}});

    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(owned, "test");
    const auto result = nlohmann::json::parse(output.str());

    CHECK(result.at("results").at(0).at("input").at("atom_ids") == nlohmann::json{"H", 9});
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
        .identity = {.source = "hydrogen.json", .record_id = "hydrogen"},
        .import_metadata = std::nullopt}};
    const auto make_result = [&records](calculation::ExecutionResult result) {
        if (result.status == calculation::ExecutionStatus::success &&
            !result.effective.has_value()) {
            result.effective = calculation::EffectiveCalculation{
                .method_id = "formal", .execution_policy = calculation::ExecutionPolicy{}};
        }
        return adapters::make_charge_calculation_result(records, {}, std::move(result));
    };

    CHECK_THROWS_AS(make_result(calculation::ExecutionResult{}), std::invalid_argument);
    CHECK_THROWS_AS(
        make_result(calculation::ExecutionResult{
            .charges = charges::ChargeSet{"formal",
                                          {{.target = {.molecule_index = 0},
                                            .charges = charges::AtomicCharges{{0.0}}}}}}),
        std::invalid_argument);
    CHECK_THROWS_AS(
        make_result(calculation::ExecutionResult{
            .charges = charges::ChargeSet{"formal",
                                          {{.target = {.molecule_index = 1},
                                            .charges = charges::AtomicCharges{{0.0, 0.0}}}}}}),
        std::invalid_argument);
    CHECK_THROWS_AS(
        make_result(calculation::ExecutionResult{
            .charges = charges::ChargeSet{"formal",
                                          {{.target = {.molecule_index = 0, .conformer_index = 1},
                                            .charges = charges::AtomicCharges{{0.0, 0.0}}}}}}),
        std::invalid_argument);
    CHECK_THROWS_AS(
        make_result(calculation::ExecutionResult{
            .charges = charges::ChargeSet{"formal",
                                          {{.target = {.molecule_index = 0},
                                            .charges = charges::AtomicCharges{{0.0, 0.0}}},
                                           {.target = {.molecule_index = 0, .conformer_index = 0},
                                            .charges = charges::AtomicCharges{{0.0, 0.0}}}}}}),
        std::invalid_argument);
    CHECK_THROWS_AS(
        make_result(calculation::ExecutionResult{
            .status = calculation::ExecutionStatus::numerical_failure,
            .charges = charges::ChargeSet{"formal",
                                          {{.target = {.molecule_index = 0},
                                            .charges = charges::AtomicCharges{{0.0, 0.0}}}}}}),
        std::invalid_argument);

    auto invalid_mapping = records;
    invalid_mapping[0].import_metadata = adapters::MoleculeImportMetadata{
        .format = adapters::MolecularSourceFormat::molecule_json,
        .atoms = {{.position = 0}},
        .conformers = {{.position = 0, .sites = {{.position = 0}}}}};
    CHECK_THROWS_AS(
        adapters::make_charge_calculation_result(
            std::move(invalid_mapping), {},
            calculation::ExecutionResult{
                .charges = charges::ChargeSet{"formal",
                                              {{.target = {.molecule_index = 0},
                                                .charges = charges::AtomicCharges{{0.0, 0.0}}}}},
                .effective =
                    calculation::EffectiveCalculation{.method_id = "formal",
                                                      .execution_policy =
                                                          calculation::ExecutionPolicy{}}}),
        std::invalid_argument);
}

TEST_CASE("result assembly requires canonical assignment order", "[adapters][json]") {
    const auto records =
        std::vector{adapters::ImportedMoleculeRecord{
                        .molecule = chargefw::core::Molecule{{chargefw::core::Atom{1}}}},
                    adapters::ImportedMoleculeRecord{
                        .molecule = chargefw::core::Molecule{{chargefw::core::Atom{1}}}}};
    const auto make_result = [&records](charges::ChargeSet charge_set) {
        return adapters::make_charge_calculation_result(
            records, {},
            {.charges = std::move(charge_set),
             .effective = calculation::EffectiveCalculation{
                 .method_id = "formal", .execution_policy = calculation::ExecutionPolicy{}}});
    };

    CHECK_NOTHROW(make_result(charges::ChargeSet{
        "formal",
        {{.target = {.molecule_index = 0}, .charges = charges::AtomicCharges{{0.0}}},
         {.target = {.molecule_index = 1}, .charges = charges::AtomicCharges{{0.0}}}}}));
    CHECK_THROWS_AS(
        make_result(charges::ChargeSet{
            "formal",
            {{.target = {.molecule_index = 1}, .charges = charges::AtomicCharges{{0.0}}},
             {.target = {.molecule_index = 0}, .charges = charges::AtomicCharges{{0.0}}}}}),
        std::invalid_argument);
}
