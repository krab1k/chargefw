#include <chargefw/adapters/charge_result.h>
#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/adapters/native/json_output.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <snitch/snitch.hpp>

#include "support/test_parameters.h"

#include <nlohmann/json.hpp>

#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace adapters = chargefw::adapters;
namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace mmcif_input = chargefw::adapters::gemmi::mmcif_input;
namespace json_output = chargefw::adapters::native::json_output;

namespace {

auto make_component_record(const std::string_view id,
                           const std::vector<std::optional<std::string>>& label_components,
                           const std::vector<std::optional<std::string>>& author_components = {})
    -> adapters::ImportedMoleculeRecord {
    auto conformers = std::vector<core::Conformer>{};
    auto metadata_conformers = std::vector<adapters::SourceConformerReference>{};
    for (std::size_t conformer_index = 0; conformer_index < label_components.size();
         ++conformer_index) {
        const auto author_component = author_components.empty()
                                          ? std::optional<std::string>{}
                                          : author_components.at(conformer_index);
        auto source = std::optional<adapters::SourceStructuralLabels>{};
        if (label_components[conformer_index].has_value() || author_component.has_value()) {
            source = adapters::SourceStructuralLabels{
                .author = {.residue = author_component},
                .label = {.residue = label_components[conformer_index]}};
        }
        auto sites = std::vector<adapters::SourceAtomReference>{
            {.position = 1, .structural_labels = std::move(source)}, {.position = 0}};
        if (conformer_index == 0) {
            metadata_conformers.push_back({.position = conformer_index, .sites = sites});
        } else {
            metadata_conformers.push_back({.position = conformer_index, .sites = std::move(sites)});
        }
        conformers.emplace_back(
            std::vector{core::Position{0.0, 0.0, static_cast<double>(conformer_index)},
                        core::Position{1.0, 0.0, static_cast<double>(conformer_index)}},
            "model-" + std::to_string(conformer_index));
    }
    auto metadata =
        adapters::MoleculeImportMetadata{.format = adapters::MolecularSourceFormat::mmcif,
                                         .atoms = metadata_conformers.front().sites,
                                         .conformers = std::move(metadata_conformers)};
    return {.molecule = core::Molecule{std::vector{core::Atom{6}, core::Atom{1}},
                                       {},
                                       std::move(conformers),
                                       std::string{id}},
            .identity = {.source = "component-fixture.cif", .record_index = 0},
            .import_metadata = std::move(metadata)};
}

auto make_component_eem_parameters() -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{
            .id = "json-eem", .method_id = "eem", .name = "JSON EEM"},
        chargefw::parameters::CommonParameters{{{.name = "kappa", .value = 2.5}}},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 5.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 9.0}}}}}};
}

} // namespace

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
    CHECK_FALSE(effective.contains("fixed_charge_groups"));

    const auto& calculated = result.at("results").at(0);
    CHECK(calculated.at("status") == "success");
    CHECK(calculated.at("input").at("record_id") == 42);
    const auto& imported = calculated.at("input").at("import");
    CHECK(imported.at("format") == "sdf");
    CHECK(imported.at("policy").at("conformer_selection") == "all");
    CHECK(imported.at("source_connectivity") == "present");
    CHECK_FALSE(imported.contains("atom_mapping"));
    CHECK_FALSE(imported.contains("conformer_mapping"));
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

TEST_CASE("JSON output projects fixed-charge group provenance", "[adapters][json]") {
    const auto result = adapters::make_charge_calculation_result(
        {{.molecule =
              chargefw::core::Molecule{
                  std::vector{chargefw::core::Atom{8}, chargefw::core::Atom{1}},
                  {},
                  {chargefw::core::Conformer{{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}, "first"},
                   chargefw::core::Conformer{{{0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}}, "second"}},
                  "water"},
          .identity = {.source = "water.json", .record_index = 0}},
         {.molecule =
              chargefw::core::Molecule{std::vector{chargefw::core::Atom{6}},
                                       {},
                                       {chargefw::core::Conformer{{{4.0, 0.0, 0.0}}, "only"}},
                                       "unaffected"},
          .identity = {.source = "other.json", .record_index = 0}}},
        {},
        {.status = calculation::ExecutionStatus::success,
         .charges = charges::ChargeSet{"eem",
                                       {{.target = {.molecule_index = 0, .conformer_index = 0},
                                         .charges = charges::AtomicCharges{{-0.25, 0.25}}},
                                        {.target = {.molecule_index = 0, .conformer_index = 1},
                                         .charges = charges::AtomicCharges{{-0.25, 0.25}}},
                                        {.target = {.molecule_index = 1, .conformer_index = 0},
                                         .charges = charges::AtomicCharges{{0.0}}}}},
         .effective = calculation::EffectiveCalculation{
             .method_id = "eem",
             .execution_policy = calculation::ExecutionPolicy{},
             .fixed_charge_groups = calculation::FixedChargeGroupsProvenance{
                 .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 0.25}},
                 .charge_provenance = "measured fixed charge",
                 .charge_totals = {{.molecule_index = 0,
                                    .original_total_charge = 0.0,
                                    .active_total_charge = 0.0},
                                   {.molecule_index = 1,
                                    .original_total_charge = 0.0,
                                    .active_total_charge = 0.0}}}}});

    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(result, "test");
    const auto document = nlohmann::json::parse(output.str());
    const auto& fixed_charge_groups =
        document.at("calculation_provenance").at("effective").at("fixed_charge_groups");
    CHECK(fixed_charge_groups.at("sources") ==
          nlohmann::json::array({{{"molecule_index", 0}, {"atom_index", 1}, {"charge", 0.25}}}));
    CHECK(fixed_charge_groups.at("charge_provenance") == "measured fixed charge");
    CHECK_FALSE(fixed_charge_groups.contains("components"));
    CHECK(
        fixed_charge_groups.at("charge_totals") ==
        nlohmann::json::array(
            {{{"molecule_index", 0}, {"original_total_charge", 0.0}, {"active_total_charge", 0.0}},
             {{"molecule_index", 1},
              {"original_total_charge", 0.0},
              {"active_total_charge", 0.0}}}));
}

TEST_CASE("JSON fixed-charge components group labeled sources without inferring unlabeled ones",
          "[adapters][json]") {
    auto records = std::vector<adapters::ImportedMoleculeRecord>{
        make_component_record("ca-first", {"CA", "CA"}),
        make_component_record("ca-second", {"CA"}),
        make_component_record("ca-other-charge", {"CA"}),
        make_component_record("author-fallback", {std::string{}}, {"MG"}),
        make_component_record("unlabeled", {std::nullopt}),
        make_component_record("conflicting-labels", {"CA", "MG"})};
    const std::vector<double> source_charges{0.4, 0.4, 0.5, 0.4, 0.4, 0.4};
    auto sources = std::vector<calculation::FixedAtomCharge>{};
    auto totals = std::vector<calculation::FixedChargeGroupChargeTotals>{};
    auto assignments = std::vector<charges::ChargeAssignment>{};
    for (std::size_t molecule_index = 0; molecule_index < records.size(); ++molecule_index) {
        const auto charge = source_charges[molecule_index];
        sources.push_back({molecule_index, 0, charge});
        totals.push_back({molecule_index, 0.0, -charge});
        for (std::size_t conformer_index = 0;
             conformer_index < records[molecule_index].molecule.conformer_count();
             ++conformer_index) {
            assignments.push_back(
                {.target = {.molecule_index = molecule_index, .conformer_index = conformer_index},
                 .charges = charges::AtomicCharges{{charge, -charge}}});
        }
    }
    const auto result = adapters::make_charge_calculation_result(
        std::move(records), {},
        {.charges = charges::ChargeSet{"eem", std::move(assignments), "component-groups"},
         .effective = calculation::EffectiveCalculation{
             .method_id = "eem",
             .parameter_set_id = "component-groups",
             .execution_policy = calculation::ExecutionPolicy{},
             .fixed_charge_groups = calculation::FixedChargeGroupsProvenance{
                 .sources = std::move(sources),
                 .charge_provenance = "explicit monatomic sources",
                 .charge_totals = std::move(totals)}}});

    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(result, "test");
    const auto document = nlohmann::json::parse(output.str());
    const auto& fixed_charge_groups =
        document.at("calculation_provenance").at("effective").at("fixed_charge_groups");
    REQUIRE(fixed_charge_groups.at("sources").size() == 6);
    const auto& components = fixed_charge_groups.at("components");
    REQUIRE(components.size() == 3);
    CHECK(components[0] == nlohmann::json{{"component_id", "CA"},
                                          {"charge_per_instance", 0.4},
                                          {"instances",
                                           {{{"molecule_index", 0}, {"atom_indices", {0}}},
                                            {{"molecule_index", 1}, {"atom_indices", {0}}}}}});
    CHECK(components[1] ==
          nlohmann::json{{"component_id", "CA"},
                         {"charge_per_instance", 0.5},
                         {"instances", {{{"molecule_index", 2}, {"atom_indices", {0}}}}}});
    CHECK(components[2] ==
          nlohmann::json{{"component_id", "MG"},
                         {"charge_per_instance", 0.4},
                         {"instances", {{{"molecule_index", 3}, {"atom_indices", {0}}}}}});
}

TEST_CASE("JSON component grouping reads labels from Gemmi imports through EEM facade results",
          "[adapters][json]") {
    const auto document = R"cif(data_first
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.label_entity_id
_atom_site.pdbx_PDB_model_num
HETATM 1 H H1 . LIG A 1 ? 0 0 0 1 20 0 1 LIG A H1 E1 1
HETATM 2 O O1 . LIG A 1 ? 2 0 0 1 20 0 1 LIG A O1 E1 1
HETATM 3 Mg MG . MG B 1 ? 0 3 0 1 20 2 1 MG B MG E2 1
#
data_second
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.label_entity_id
_atom_site.pdbx_PDB_model_num
HETATM 1 H H1 . LIG A 1 ? 0 0 1 1 20 0 1 LIG A H1 E1 1
HETATM 2 O O1 . LIG A 1 ? 2 0 1 1 20 0 1 LIG A O1 E1 1
HETATM 3 Mg MG . MG B 1 ? 0 3 1 1 20 2 1 MG B MG E2 1
#
)cif";
    auto input = std::istringstream{document};
    auto reader = mmcif_input::MmcifReader{input, "component-input.cif"};
    auto records = std::vector<adapters::ImportedMoleculeRecord>{};
    while (auto record = reader.next()) {
        records.push_back(std::move(*record));
    }
    REQUIRE(records.size() == 2);
    REQUIRE(records[0].import_metadata.has_value());
    REQUIRE(records[0].import_metadata->atoms[2].structural_labels.has_value());
    CHECK(records[0].import_metadata->atoms[2].structural_labels->label.residue == "MG");

    auto molecules = std::vector<core::Molecule>{};
    for (const auto& record : records) {
        molecules.push_back(record.molecule);
    }
    const auto parameters = make_component_eem_parameters();
    auto execution = calculation::calculate(calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::move(molecules)},
        .parameter_sets = {parameters},
        .method_id = "eem",
        .fixed_charge_groups = calculation::FixedChargeGroups{
            .sources = {{.molecule_index = 0, .atom_index = 2, .charge = 0.4},
                        {.molecule_index = 1, .atom_index = 2, .charge = 0.4}},
            .charge_provenance = "imported fixed ions"}}));
    REQUIRE(execution.calculated());
    const auto result = adapters::make_charge_calculation_result(
        std::move(records), {.method_id = "eem", .parameter_set_id = "json-eem"},
        std::move(execution));

    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(result, "test");
    const auto json = nlohmann::json::parse(output.str());
    const auto& fixed_charge_groups =
        json.at("calculation_provenance").at("effective").at("fixed_charge_groups");
    const auto& charge_totals = fixed_charge_groups.at("charge_totals");
    REQUIRE(charge_totals.size() == 2);
    for (std::size_t molecule_index = 0; molecule_index < charge_totals.size(); ++molecule_index) {
        CHECK(charge_totals[molecule_index].at("original_total_charge") == 2.0);
        CHECK(charge_totals[molecule_index].at("active_total_charge") == 0.0);
        auto prescribed_total = 0.0;
        for (const auto& source : fixed_charge_groups.at("sources")) {
            if (source.at("molecule_index") == molecule_index) {
                prescribed_total += source.at("charge").get<double>();
            }
        }
        CHECK(charge_totals[molecule_index].at("active_total_charge").get<double>() +
                  prescribed_total ==
              0.4);
    }
    REQUIRE(fixed_charge_groups.at("components").size() == 1);
    CHECK(fixed_charge_groups.at("components")[0].at("component_id") == "MG");
    CHECK(fixed_charge_groups.at("components")[0].at("charge_per_instance") == 0.4);
    REQUIRE(fixed_charge_groups.at("components")[0].at("instances").size() == 2);
    CHECK(fixed_charge_groups.at("components")[0].at("instances")[0] ==
          nlohmann::json{{"molecule_index", 0}, {"atom_indices", {2}}});
    CHECK(fixed_charge_groups.at("components")[0].at("instances")[1] ==
          nlohmann::json{{"molecule_index", 1}, {"atom_indices", {2}}});
}

TEST_CASE("JSON output serializes a cancelled result without assignments", "[adapters][json]") {
    const auto owned = adapters::make_charge_calculation_result(
        {{.molecule = chargefw::core::Molecule{std::vector{chargefw::core::Atom{8},
                                                           chargefw::core::Atom{1}}},
          .identity = {.source = "water.sdf", .record_index = 0}}},
        {},
        calculation::ExecutionResult{
            .status = calculation::ExecutionStatus::cancelled,
            .effective = calculation::EffectiveCalculation{
                .method_id = "eem",
                .execution_policy = calculation::ExecutionPolicy{},
                .fixed_charge_groups = calculation::FixedChargeGroupsProvenance{
                    .sources = {{.molecule_index = 0, .atom_index = 0, .charge = 0.5}},
                    .charge_provenance = "",
                    .charge_totals = {{.molecule_index = 0,
                                       .original_total_charge = 0.0,
                                       .active_total_charge = 0.0}}}}});

    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(owned, "test");
    const auto result = nlohmann::json::parse(output.str());

    CHECK(result.at("status") == "cancelled");
    CHECK(result.at("calculation_provenance").at("requested").at("execution").at("kind") == "auto");
    CHECK(result.at("diagnostics").at(0).at("code") == "calculation_cancelled");
    const auto& record = result.at("results").at(0);
    CHECK(record.at("status") == "cancelled");
    CHECK_FALSE(record.contains("assignments"));
    CHECK(result.at("calculation_provenance")
              .at("effective")
              .at("fixed_charge_groups")
              .at("charge_provenance")
              .get<std::string>()
              .empty());
}

TEST_CASE("JSON output retains fixed-charge group metadata on numerical failure without charges",
          "[adapters][json]") {
    const auto result = adapters::make_charge_calculation_result(
        {{.molecule = chargefw::core::Molecule{std::vector{chargefw::core::Atom{6},
                                                           chargefw::core::Atom{1}}},
          .identity = {.source = "carbon.json", .record_index = 0}}},
        {},
        calculation::ExecutionResult{
            .status = calculation::ExecutionStatus::numerical_failure,
            .effective = calculation::EffectiveCalculation{
                .method_id = "eem",
                .execution_policy = calculation::ExecutionPolicy{},
                .fixed_charge_groups = calculation::FixedChargeGroupsProvenance{
                    .sources = {{.molecule_index = 0, .atom_index = 0, .charge = -0.2}},
                    .charge_provenance = "fixed value",
                    .charge_totals = {{.molecule_index = 0,
                                       .original_total_charge = 0.0,
                                       .active_total_charge = 0.0}}}}});
    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(result, "test");
    const auto document = nlohmann::json::parse(output.str());

    CHECK(document.at("status") == "numerical_failure");
    CHECK(document.at("results").at(0).at("status") == "numerical_failure");
    CHECK_FALSE(document.at("results").at(0).contains("assignments"));
    CHECK(document.at("calculation_provenance")
              .at("effective")
              .at("fixed_charge_groups")
              .contains("charge_totals"));
}

TEST_CASE("JSON output retains cancelled results without effective provenance",
          "[adapters][json]") {
    const auto result = adapters::make_charge_calculation_result(
        {{.molecule = chargefw::core::Molecule{std::vector{chargefw::core::Atom{8}}},
          .identity = {.source = "water.sdf", .record_index = 0}}},
        {}, calculation::ExecutionResult{.status = calculation::ExecutionStatus::cancelled});
    auto output = std::ostringstream{};
    json_output::JsonWriter{output}.write(result, "test");
    const auto document = nlohmann::json::parse(output.str());

    CHECK(document.at("status") == "cancelled");
    CHECK_FALSE(
        document.at("calculation_provenance").at("effective").contains("fixed_charge_groups"));
}

TEST_CASE("result assembly validates fixed-charge group provenance", "[adapters][json]") {
    const auto records = std::vector{
        adapters::ImportedMoleculeRecord{.molecule = chargefw::core::Molecule{std::vector{
                                             chargefw::core::Atom{6}, chargefw::core::Atom{8}}}}};
    const auto make_result =
        [&records](calculation::FixedChargeGroupsProvenance fixed_charge_groups) {
            return adapters::make_charge_calculation_result(
                records, {},
                {.status = calculation::ExecutionStatus::cancelled,
                 .effective = calculation::EffectiveCalculation{
                     .method_id = "eem",
                     .execution_policy = calculation::ExecutionPolicy{},
                     .fixed_charge_groups = std::move(fixed_charge_groups)}});
        };
    const auto valid = calculation::FixedChargeGroupsProvenance{
        .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 0.5}},
        .charge_provenance = "caller label",
        .charge_totals = {
            {.molecule_index = 0, .original_total_charge = 0.0, .active_total_charge = 0.0}}};
    CHECK_NOTHROW(make_result(valid));

    auto invalid = valid;
    invalid.sources.clear();
    CHECK_THROWS_AS(make_result(invalid), std::invalid_argument);
    invalid = valid;
    invalid.sources[0].molecule_index = 1;
    CHECK_THROWS_AS(make_result(invalid), std::invalid_argument);
    invalid = valid;
    invalid.sources[0].atom_index = 2;
    CHECK_THROWS_AS(make_result(invalid), std::invalid_argument);
    invalid = valid;
    invalid.sources[0].charge = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(make_result(invalid), std::invalid_argument);
    invalid = valid;
    invalid.sources.push_back(invalid.sources.front());
    CHECK_THROWS_AS(make_result(invalid), std::invalid_argument);
    invalid = valid;
    invalid.charge_totals.clear();
    CHECK_THROWS_AS(make_result(invalid), std::invalid_argument);
    invalid = valid;
    invalid.charge_totals[0].molecule_index = 1;
    CHECK_THROWS_AS(make_result(invalid), std::invalid_argument);
    invalid = valid;
    invalid.charge_totals[0].active_total_charge = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(make_result(invalid), std::invalid_argument);
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
    CHECK(result.at("calculation_provenance").at("requested").at("execution").at("kind") == "auto");
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
