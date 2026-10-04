#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <snitch/snitch.hpp>

#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace calculation = chargefw::calculation;
namespace core = chargefw::core;
namespace methods = chargefw::methods;

namespace {

auto make_isolated_ion_pair() -> core::Molecule {
    return core::Molecule{
        std::vector{core::Atom{1}, core::Atom{12, 2}},
        {},
        {core::Conformer{{core::Position{0.0, 0.0, 0.0}, core::Position{3.0, 0.0, 0.0}}}},
        "isolated-ion-pair"};
}

auto make_interleaved_hydrogen_magnesium_oxygen(
    const core::Position oxygen_position = core::Position{1.0, 0.0, 0.0}) -> core::Molecule {
    return core::Molecule{std::vector{core::Atom{1, 0, "active-H"}, core::Atom{12, 2, "fixed-Mg"},
                                      core::Atom{8, 0, "active-O"}},
                          {},
                          {core::Conformer{{core::Position{0.0, 0.0, 0.0},
                                            core::Position{3.0, 0.0, 0.0}, oxygen_position}}},
                          "interleaved-H-Mg-O"};
}

auto make_bonded_pair() -> core::Molecule {
    return core::Molecule{std::vector{core::Atom{1}, core::Atom{6}},
                          std::vector{core::Bond{0, 1}},
                          {},
                          "bonded-pair"};
}

auto make_molecule_with_positions(std::vector<core::Atom> atoms,
                                  std::vector<std::vector<core::Position>> conformers,
                                  std::vector<core::Bond> bonds = {}) -> core::Molecule {
    auto conformer_values = std::vector<core::Conformer>{};
    conformer_values.reserve(conformers.size());
    for (auto& positions : conformers) {
        conformer_values.emplace_back(std::move(positions));
    }
    return core::Molecule{std::move(atoms), std::move(bonds), std::move(conformer_values)};
}

auto assessment_error(calculation::AssessmentRequest request) -> std::string {
    try {
        static_cast<void>(calculation::assess(std::move(request)));
    } catch (const std::invalid_argument& error) {
        return error.what();
    }
    return {};
}

auto make_eem_parameters() -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{
            .id = "test-eem", .method_id = "eem", .name = "Test EEM"},
        chargefw::parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 10.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 10.0}}}}}};
}

auto make_embedding_eem_parameters() -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{
            .id = "embedding-eem", .method_id = "eem", .name = "Embedding EEM"},
        chargefw::parameters::CommonParameters{{{.name = "kappa", .value = 2.5}}},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 5.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 9.0}}}}}};
}

auto make_two_conformer_interleaved_pair(const int hydrogen_charge = 0,
                                         const int magnesium_charge = 2) -> core::Molecule {
    return core::Molecule{
        std::vector{core::Atom{1, hydrogen_charge}, core::Atom{12, magnesium_charge},
                    core::Atom{8, 0}},
        {},
        {core::Conformer{{{0.0, 0.0, 0.0}, {0.0, 3.0, 0.0}, {2.0, 0.0, 0.0}}, "first"},
         core::Conformer{{{0.0, 0.0, 0.0}, {4.0, 1.0, 0.0}, {2.0, 0.0, 0.0}}, "second"}},
        "interleaved-pair"};
}

auto make_sqeqp_embedding_parameters() -> chargefw::parameters::ParameterSet {
    auto atoms = std::vector<chargefw::parameters::AtomParameterEntry>{
        {.key = chargefw::test::plain_atom_key(1),
         .parameters = {{.name = "electronegativity", .value = 4.528},
                        {.name = "hardness", .value = 13.8904},
                        {.name = "width", .value = 1.0},
                        {.name = "q0", .value = 0.3}}},
        {.key = chargefw::test::plain_atom_key(8),
         .parameters = {{.name = "electronegativity", .value = 8.741},
                        {.name = "hardness", .value = 13.364},
                        {.name = "width", .value = -0.5},
                        {.name = "q0", .value = 0.1}}}};
    auto bonds = std::vector<chargefw::parameters::BondParameterEntry>{
        {.key = chargefw::test::single_bond_key(1, 8),
         .parameters = {{.name = "kappa", .value = 1.0}}}};
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{
            .id = "embedding-sqeqp", .method_id = "sqeqp", .name = "Embedding SQE+qp"},
        {},
        chargefw::parameters::AtomParameters{std::move(atoms)},
        chargefw::parameters::BondParameters{std::move(bonds)}};
}

auto make_sqeqp_embedding_molecule(const int magnesium_charge = 1) -> core::Molecule {
    return core::Molecule{
        std::vector{core::Atom{1, -1}, core::Atom{12, magnesium_charge}, core::Atom{8, 0}},
        {core::Bond{0, 2}},
        {core::Conformer{{{0.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {1.5, 0.0, 0.0}}, "ion-interleaved"}},
        "sqeqp-embedding"};
}

auto make_sqe_active_molecule(const int hydrogen_charge, const int magnesium_charge)
    -> core::Molecule {
    return core::Molecule{
        std::vector{core::Atom{1, hydrogen_charge}, core::Atom{12, magnesium_charge},
                    core::Atom{8, 0}},
        {core::Bond{0, 2}},
        {core::Conformer{{{0.0, 0.0, 0.0}, {0.0, 3.0, 0.0}, {1.5, 0.0, 0.0}}, "active"}},
        "sqe-active"};
}

auto make_fractional_sqe_source_molecule() -> core::Molecule {
    return core::Molecule{
        std::vector{core::Atom{1, 0}, core::Atom{12, 1}, core::Atom{8, 0}, core::Atom{11, -1}},
        {core::Bond{0, 2}},
        {core::Conformer{{{0.0, 0.0, 0.0}, {0.0, 3.0, 0.0}, {1.5, 0.0, 0.0}, {4.0, 2.0, 0.0}},
                         "fractional-sources"}},
        "fractional-sqe-sources"};
}

auto make_sqe_embedding_parameters(const std::string_view method_id)
    -> chargefw::parameters::ParameterSet {
    auto atoms = std::vector<chargefw::parameters::AtomParameterEntry>{
        {.key = chargefw::test::plain_atom_key(1),
         .parameters = {{.name = "electronegativity", .value = 4.528},
                        {.name = "hardness", .value = 13.8904},
                        {.name = "width", .value = 1.0}}},
        {.key = chargefw::test::plain_atom_key(8),
         .parameters = {{.name = "electronegativity", .value = 8.741},
                        {.name = "hardness", .value = 13.364},
                        {.name = "width", .value = 1.0}}}};
    auto bonds = std::vector<chargefw::parameters::BondParameterEntry>{
        {.key = chargefw::test::single_bond_key(1, 8),
         .parameters = {{.name = "kappa", .value = 1.0}}}};
    const auto method = std::string{method_id};
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{
            .id = "embedding-" + method, .method_id = method, .name = "Embedding " + method},
        {},
        chargefw::parameters::AtomParameters{std::move(atoms)},
        chargefw::parameters::BondParameters{std::move(bonds)}};
}

auto make_hydrogen_only_eem_parameters() -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{
            .id = "test-eem-h-only", .method_id = "eem", .name = "Test EEM H-only"},
        chargefw::parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 10.0}}}}}};
}

} // namespace

TEST_CASE("explicit unsupported execution has no selected plan or fallback",
          "[calculation][planning]") {
    auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "formal",
        .execution_selection = calculation::ExecutionSelection{
            calculation::ExecutionSelectionKind::cover, calculation::minimum_reduced_radius}});

    CHECK(assessment.plans().empty());
    REQUIRE(assessment.rejections().size() == 1);
    REQUIRE(assessment.rejections()[0].policy.has_value());
    CHECK(assessment.rejections()[0].policy->mode() == calculation::ExecutionMode::cover);
    REQUIRE(assessment.rejections()[0].issues.size() == 1);
    CHECK(std::get<methods::ExecutionIssue>(assessment.rejections()[0].issues[0]).kind ==
          methods::ExecutionIssueKind::unsupported_execution_mode);

    const auto result = calculation::calculate(assessment);
    CHECK(result.status == calculation::ExecutionStatus::no_executable_plan);
    CHECK_FALSE(result.calculated());
}

TEST_CASE("explicit no-plan assessment reports rejected scientific prerequisites",
          "[calculation][planning]") {
    auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "smpqeq"});

    CHECK(assessment.plans().empty());
    REQUIRE(assessment.rejections().size() == 1);
    CHECK_FALSE(assessment.rejections()[0].policy.has_value());
    CHECK(assessment.rejections()[0].method_id == "smpqeq");
    CHECK_FALSE(assessment.rejections()[0].issues.empty());

    const auto result = calculation::calculate(assessment);
    CHECK(result.status == calculation::ExecutionStatus::no_executable_plan);
    CHECK_FALSE(result.calculated());
}

TEST_CASE("valid fixed charge selectors reach EEM planning with active molecules",
          "[calculation][planning]") {
    auto request = calculation::AssessmentRequest{
        .molecules =
            core::MoleculeCollection{std::vector{make_interleaved_hydrogen_magnesium_oxygen(),
                                                 make_interleaved_hydrogen_magnesium_oxygen(),
                                                 make_interleaved_hydrogen_magnesium_oxygen()}},
        .parameter_sets = {make_eem_parameters()},
        .method_id = "eem",
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 2, .atom_index = 1, .charge = -0.5},
                        {.molecule_index = 1, .atom_index = 1, .charge = 0.0},
                        {.molecule_index = 0, .atom_index = 1, .charge = 0.25}},
            .charge_provenance = "test fixed ions"}};
    const auto assessment = calculation::assess(std::move(request));

    REQUIRE(assessment.plans().size() == 1);
    CHECK(assessment.default_plan()->candidate().method->id() == "eem");
    CHECK(assessment.default_plan()->policy().mode() == calculation::ExecutionMode::full);
    REQUIRE(assessment.molecules().size() == 3);
    CHECK(assessment.molecules()[1].atom(1).atomic_number() == 12);
    CHECK(assessment.molecules()[1].atom(1).name() == "fixed-Mg");
    REQUIRE(assessment.rejections().size() == 2);
    for (const auto& rejection : assessment.rejections()) {
        CHECK(rejection.method_id == "eem");
        REQUIRE(rejection.policy.has_value());
        REQUIRE(rejection.issues.size() == 1);
        CHECK(std::get<methods::ExecutionIssue>(rejection.issues[0]).kind ==
              methods::ExecutionIssueKind::unsupported_execution_mode);
    }
    const auto execution = calculation::calculate(assessment, 2);
    REQUIRE(execution.calculated());
    REQUIRE(execution.charges->size() == 3);
    REQUIRE(execution.effective.has_value());
    REQUIRE(execution.effective->fixed_charge_embedding.has_value());
    const auto& provenance = *execution.effective->fixed_charge_embedding;
    CHECK(provenance.charge_provenance == "test fixed ions");
    REQUIRE(provenance.sources.size() == 3);
    CHECK(provenance.sources[0].molecule_index == 0);
    CHECK(provenance.sources[0].atom_index == 1);
    CHECK(provenance.sources[0].charge == 0.25);
    CHECK(provenance.sources[1].charge == 0.0);
    CHECK(provenance.sources[2].charge == -0.5);
    REQUIRE(provenance.charge_totals.size() == 3);
    CHECK(provenance.charge_totals[0].original_total_charge == 2.0);
    CHECK(provenance.charge_totals[0].active_total_charge == 0.0);
    CHECK(provenance.charge_totals[1].original_total_charge == 2.0);
    CHECK(provenance.charge_totals[1].active_total_charge == 0.0);
    CHECK(provenance.charge_totals[2].original_total_charge == 2.0);
    CHECK(provenance.charge_totals[2].active_total_charge == 0.0);
    CHECK(execution.charges->assignment(0).charges[1] == 0.25);
    CHECK(execution.charges->assignment(1).charges[1] == 0.0);
    CHECK(execution.charges->assignment(2).charges[1] == -0.5);
}

TEST_CASE("native embedded EEM executes, reassembles, and retains provenance",
          "[calculation][planning]") {
    auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{
            make_two_conformer_interleaved_pair(),
            core::Molecule{std::vector{core::Atom{1, -1}},
                           {},
                           {core::Conformer{{{8.0, 0.0, 0.0}}, "unaffected"}},
                           "unaffected"}}},
        .parameter_sets = {make_embedding_eem_parameters()},
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 0.4}},
            .charge_provenance = "fractional fixed magnesium"}});
    auto moved_assessment = std::move(assessment);

    REQUIRE(moved_assessment.plans().size() == 1);
    const auto& plan = moved_assessment.plans()[0];
    CHECK(plan.candidate().method->id() == "eem");
    CHECK(plan.policy().mode() == calculation::ExecutionMode::full);

    const auto serial = calculation::calculate(moved_assessment, plan, 1);
    const auto parallel = calculation::calculate(moved_assessment, plan, 2);
    const auto repeated = calculation::calculate(moved_assessment, plan, 1);
    for (const auto* result : {&serial, &parallel, &repeated}) {
        REQUIRE(result->calculated());
        REQUIRE(result->charges->size() == 3);
        REQUIRE(result->effective.has_value());
        REQUIRE(result->effective->fixed_charge_embedding.has_value());
    }

    const auto expected_pair = [](const double source_potential_h,
                                  const double source_potential_o) {
        constexpr auto active_charge = 0.0;
        constexpr auto cross_interaction = 2.5 / 2.0;
        const auto rhs_difference = -1.0 - source_potential_h + 2.0 + source_potential_o;
        const auto hydrogen_charge = (rhs_difference + (9.0 - cross_interaction) * active_charge) /
                                     (5.0 + 9.0 - 2.0 * cross_interaction);
        return std::array{hydrogen_charge, active_charge - hydrogen_charge};
    };
    const auto first_expected = expected_pair(2.5 * 0.4 / 3.0, 2.5 * 0.4 / std::sqrt(13.0));
    const auto second_expected =
        expected_pair(2.5 * 0.4 / std::sqrt(17.0), 2.5 * 0.4 / std::sqrt(5.0));
    const auto no_field = expected_pair(0.0, 0.0);
    for (std::size_t index = 0; index < serial.charges->size(); ++index) {
        const auto& assignment = serial.charges->assignment(index);
        CHECK(std::ranges::equal(assignment.charges.values(),
                                 parallel.charges->assignment(index).charges.values()));
        CHECK(std::ranges::equal(assignment.charges.values(),
                                 repeated.charges->assignment(index).charges.values()));
    }
    CHECK(serial.charges->assignment(0).target.molecule_index == 0);
    CHECK(serial.charges->assignment(0).target.conformer_index == 0);
    CHECK(serial.charges->assignment(1).target.molecule_index == 0);
    CHECK(serial.charges->assignment(1).target.conformer_index == 1);
    CHECK(serial.charges->assignment(2).target.molecule_index == 1);
    CHECK(serial.charges->assignment(2).target.conformer_index == 0);

    for (std::size_t conformer = 0; conformer < 2; ++conformer) {
        const auto& expected = conformer == 0 ? first_expected : second_expected;
        const auto& values = serial.charges->assignment(conformer).charges;
        CHECK(std::abs(values[0] - expected[0]) < 1e-10);
        CHECK(std::abs(values[0] - no_field[0]) > 1e-4);
        CHECK(values[1] == 0.4);
        CHECK(std::abs(values[2] - expected[1]) < 1e-10);
        CHECK(std::abs(values.total() - 0.4) < 1e-12);
        CHECK(std::abs(values[0] + values[2]) < 1e-12);
    }
    CHECK(serial.charges->assignment(2).charges[0] == -1.0);
    CHECK(serial.effective->fixed_charge_embedding->charge_provenance ==
          "fractional fixed magnesium");
    REQUIRE(serial.effective->fixed_charge_embedding->sources.size() == 1);
    CHECK(serial.effective->fixed_charge_embedding->sources[0].molecule_index == 0);
    CHECK(serial.effective->fixed_charge_embedding->sources[0].atom_index == 1);
    CHECK(serial.effective->fixed_charge_embedding->sources[0].charge == 0.4);
    REQUIRE(serial.effective->fixed_charge_embedding->charge_totals.size() == 2);
    CHECK(serial.effective->fixed_charge_embedding->charge_totals[0].original_total_charge == 2.0);
    CHECK(serial.effective->fixed_charge_embedding->charge_totals[0].active_total_charge == 0.0);
    CHECK(serial.effective->fixed_charge_embedding->charge_totals[1].original_total_charge == -1.0);
    CHECK(serial.effective->fixed_charge_embedding->charge_totals[1].active_total_charge == -1.0);
    CHECK(moved_assessment.molecules()[0].atom(1).formal_charge() == 2);

    auto warning_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_two_conformer_interleaved_pair()}},
        .parameter_sets = {make_embedding_eem_parameters()},
        .method_id = "eem",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
        .resource_policy = {.cutoff_atom_threshold = 1},
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 0.4}}}});
    REQUIRE(warning_assessment.default_plan() != nullptr);
    REQUIRE(warning_assessment.default_plan()->warnings().size() == 1);
    CHECK(warning_assessment.default_plan()->warnings()[0].kind ==
          methods::ExecutionIssueKind::resource_threshold_exceeded);
    const auto warning_result = calculation::calculate(warning_assessment);
    REQUIRE(warning_result.calculated());
    CHECK(warning_result.effective->fixed_charge_embedding.has_value());

    auto automatic_limited_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_two_conformer_interleaved_pair()}},
        .parameter_sets = {make_embedding_eem_parameters()},
        .method_id = "eem",
        .resource_policy = {.cutoff_atom_threshold = 1},
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 0.4}}}});
    CHECK(automatic_limited_assessment.plans().empty());
    CHECK(std::ranges::any_of(automatic_limited_assessment.rejections(), [](const auto& rejection) {
        return rejection.policy.has_value() &&
               rejection.policy->mode() == calculation::ExecutionMode::full &&
               std::ranges::any_of(rejection.issues, [](const auto& issue) {
                   const auto* execution_issue = std::get_if<methods::ExecutionIssue>(&issue);
                   return execution_issue != nullptr &&
                          execution_issue->kind ==
                              methods::ExecutionIssueKind::resource_threshold_exceeded;
               });
    }));
}

TEST_CASE("native SQE+qp normalizes reference charges to prepared active formal charge",
          "[calculation][planning]") {
    auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_sqeqp_embedding_molecule()}},
        .parameter_sets = {make_sqeqp_embedding_parameters()},
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 0.3}},
            .charge_provenance = "fixed magnesium"}});
    REQUIRE(assessment.plans().size() == 1);
    const auto& plan = assessment.plans()[0];
    CHECK(plan.candidate().method->id() == "sqeqp");
    CHECK(plan.policy().mode() == calculation::ExecutionMode::full);
    const auto interaction = [](const double distance, const double width_a, const double width_b) {
        const auto width_sum = 2.0 * width_a * width_a + 2.0 * width_b * width_b;
        return width_sum == 0.0 ? 1.0 / distance
                                : std::erf(distance / std::sqrt(width_sum)) / distance;
    };
    constexpr auto seed_h = -0.4;
    constexpr auto seed_o = -0.6;
    const auto active_interaction = interaction(1.5, 1.0, -0.5);
    const auto source_rhs_h =
        -4.528 - 0.3 * interaction(2.0, 1.0, 0.0) - active_interaction * seed_o;
    const auto source_rhs_o =
        -8.741 - 0.3 * interaction(2.5, -0.5, 0.0) - active_interaction * seed_h;
    const auto no_field_rhs_h = -4.528 - active_interaction * seed_o;
    const auto no_field_rhs_o = -8.741 - active_interaction * seed_h;
    const auto denominator = 13.8904 + 13.364 - 2.0 * active_interaction + 1.0;
    const auto expected_h = seed_h + (source_rhs_h - source_rhs_o) / denominator;
    const auto expected_o = seed_o - (source_rhs_h - source_rhs_o) / denominator;
    const auto no_field_h = seed_h + (no_field_rhs_h - no_field_rhs_o) / denominator;
    const auto no_field_o = seed_o - (no_field_rhs_h - no_field_rhs_o) / denominator;

    const auto serial = calculation::calculate(assessment, plan, 1);
    const auto parallel = calculation::calculate(assessment, plan, 2);
    const auto repeated = calculation::calculate(assessment, plan, 1);
    for (const auto* result : {&serial, &parallel, &repeated}) {
        REQUIRE(result->calculated());
        REQUIRE(result->charges->size() == 1);
        REQUIRE(result->effective->fixed_charge_embedding.has_value());
        const auto& values = result->charges->assignment(0).charges;
        CHECK(values[1] == 0.3);
        CHECK(std::abs(values[0] + values[2] + 1.0) < 1e-12);
        CHECK(std::abs(values.total() + 0.7) < 1e-12);
    }
    CHECK(std::abs(serial.charges->assignment(0).charges[0] - expected_h) < 1e-12);
    CHECK(std::abs(serial.charges->assignment(0).charges[2] - expected_o) < 1e-12);
    CHECK(std::abs(serial.charges->assignment(0).charges[0] - no_field_h) > 1e-4);
    CHECK(std::abs(serial.charges->assignment(0).charges[2] - no_field_o) > 1e-4);
    CHECK(std::ranges::equal(serial.charges->assignment(0).charges.values(),
                             parallel.charges->assignment(0).charges.values()));
    CHECK(std::ranges::equal(serial.charges->assignment(0).charges.values(),
                             repeated.charges->assignment(0).charges.values()));
    const auto& provenance = *serial.effective->fixed_charge_embedding;
    REQUIRE(provenance.sources.size() == 1);
    CHECK(provenance.sources[0].atom_index == 1);
    CHECK(provenance.sources[0].charge == 0.3);
    REQUIRE(provenance.charge_totals.size() == 1);
    CHECK(provenance.charge_totals[0].original_total_charge == 0.0);
    CHECK(provenance.charge_totals[0].active_total_charge == -1.0);
    CHECK(assessment.molecules()[0].atom(1).formal_charge() == 1);

    auto cutoff_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_sqeqp_embedding_molecule()}},
        .parameter_sets = {make_sqeqp_embedding_parameters()},
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff,
                                            calculation::minimum_reduced_radius},
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 0.3}}}});
    CHECK(cutoff_assessment.plans().empty());
    CHECK(std::ranges::any_of(cutoff_assessment.rejections(), [](const auto& rejection) {
        return rejection.method_id == "sqeqp" && rejection.policy.has_value() &&
               rejection.policy->mode() == calculation::ExecutionMode::cutoff;
    }));
}

TEST_CASE("fixed source imported charge does not change the active EEM charge",
          "[calculation][planning][embedding]") {
    auto results = std::vector<calculation::ExecutionResult>{};
    for (const auto magnesium_charge : {0, 2}) {
        auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{
                make_two_conformer_interleaved_pair(0, magnesium_charge)}},
            .parameter_sets = {make_embedding_eem_parameters()},
            .method_id = "eem",
            .fixed_charge_embedding = calculation::FixedChargeEmbedding{
                .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}}}});
        const auto result = calculation::calculate(assessment);
        REQUIRE(result.calculated());
        REQUIRE(result.effective->fixed_charge_embedding.has_value());
        const auto& totals = result.effective->fixed_charge_embedding->charge_totals[0];
        CHECK(totals.original_total_charge == magnesium_charge);
        CHECK(totals.active_total_charge == 0.0);
        for (const auto& assignment : result.charges->assignments()) {
            CHECK(assignment.charges[1] == 2.0);
            CHECK(std::abs(assignment.charges[0] + assignment.charges[2]) < 1.0e-12);
            CHECK(std::abs(assignment.charges.total() - 2.0) < 1.0e-12);
        }
        results.push_back(result);
    }

    REQUIRE(results[0].charges->size() == results[1].charges->size());
    for (std::size_t index = 0; index < results[0].charges->size(); ++index) {
        CHECK(std::ranges::equal(results[0].charges->assignment(index).charges.values(),
                                 results[1].charges->assignment(index).charges.values()));
    }

    auto changed_field = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_two_conformer_interleaved_pair()}},
        .parameter_sets = {make_embedding_eem_parameters()},
        .method_id = "eem",
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 1.0}}}});
    const auto changed = calculation::calculate(changed_field);
    REQUIRE(changed.calculated());
    CHECK(changed.effective->fixed_charge_embedding->charge_totals[0].active_total_charge == 0.0);
    CHECK(std::abs(changed.charges->assignment(0).charges.total() - 1.0) < 1.0e-12);
    CHECK(std::abs(changed.charges->assignment(0).charges[0] -
                   results[1].charges->assignment(0).charges[0]) > 1.0e-4);
}

TEST_CASE("prepared active formal charge sets EEM and SQE+qp totals",
          "[calculation][planning][embedding]") {
    auto sqeqp_results = std::vector<calculation::ExecutionResult>{};
    for (const auto magnesium_charge : {0, 2}) {
        auto eem_assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{
                make_two_conformer_interleaved_pair(-1, magnesium_charge)}},
            .parameter_sets = {make_embedding_eem_parameters()},
            .method_id = "eem",
            .fixed_charge_embedding = calculation::FixedChargeEmbedding{
                .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}}}});
        const auto eem = calculation::calculate(eem_assessment);
        REQUIRE(eem.calculated());
        const auto& eem_totals = eem.effective->fixed_charge_embedding->charge_totals[0];
        CHECK(eem_totals.original_total_charge == magnesium_charge - 1);
        CHECK(eem_totals.active_total_charge == -1.0);
        for (const auto& assignment : eem.charges->assignments()) {
            CHECK(std::abs(assignment.charges[0] + assignment.charges[2] + 1.0) < 1.0e-12);
            CHECK(std::abs(assignment.charges.total() - 1.0) < 1.0e-12);
        }

        auto sqeqp_assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{
                make_sqeqp_embedding_molecule(magnesium_charge)}},
            .parameter_sets = {make_sqeqp_embedding_parameters()},
            .method_id = "sqeqp",
            .fixed_charge_embedding = calculation::FixedChargeEmbedding{
                .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}}}});
        const auto sqeqp = calculation::calculate(sqeqp_assessment);
        REQUIRE(sqeqp.calculated());
        const auto& sqeqp_totals = sqeqp.effective->fixed_charge_embedding->charge_totals[0];
        CHECK(sqeqp_totals.original_total_charge == magnesium_charge - 1);
        CHECK(sqeqp_totals.active_total_charge == -1.0);
        CHECK(std::abs(sqeqp.charges->assignment(0).charges[0] +
                       sqeqp.charges->assignment(0).charges[2] + 1.0) < 1.0e-12);
        CHECK(std::abs(sqeqp.charges->assignment(0).charges.total() - 1.0) < 1.0e-12);
        sqeqp_results.push_back(sqeqp);
    }
    CHECK(std::ranges::equal(sqeqp_results[0].charges->assignment(0).charges.values(),
                             sqeqp_results[1].charges->assignment(0).charges.values()));
}

TEST_CASE("SQE+q0 retains prepared component charges while source atoms are fixed",
          "[calculation][planning][sqeq0]") {
    auto results = std::vector<calculation::ExecutionResult>{};
    for (const auto magnesium_charge : {0, 2}) {
        auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{
                make_sqe_active_molecule(-1, magnesium_charge)}},
            .parameter_sets = {make_sqe_embedding_parameters("sqeq0")},
            .method_id = "sqeq0",
            .fixed_charge_embedding = calculation::FixedChargeEmbedding{
                .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}}}});
        const auto result = calculation::calculate(assessment);
        REQUIRE(result.calculated());
        REQUIRE(result.effective->fixed_charge_embedding.has_value());
        const auto& values = result.charges->assignment(0).charges;
        CHECK(values[1] == 2.0);
        CHECK(std::abs(values[0] + values[2] + 1.0) < 1.0e-12);
        CHECK(std::abs(values.total() - 1.0) < 1.0e-12);
        const auto& totals = result.effective->fixed_charge_embedding->charge_totals[0];
        CHECK(totals.original_total_charge == magnesium_charge - 1);
        CHECK(totals.active_total_charge == -1.0);
        results.push_back(result);
    }
    CHECK(std::ranges::equal(results[0].charges->assignment(0).charges.values(),
                             results[1].charges->assignment(0).charges.values()));

    const auto sqe_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_sqe_active_molecule(-1, 2)}},
        .method_id = "sqe",
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}}}});
    CHECK(sqe_assessment.plans().empty());
    REQUIRE(sqe_assessment.rejections().size() == 1);
    CHECK(sqe_assessment.rejections()[0].parameter_set_id == std::nullopt);
    REQUIRE(sqe_assessment.rejections()[0].issues.size() == 1);
    CHECK(std::get<methods::PrerequisiteIssue>(sqe_assessment.rejections()[0].issues[0]).kind ==
          methods::PrerequisiteIssueKind::unsupported_molecule);
}

TEST_CASE("SQE family preserves active totals with fractional fixed-source values",
          "[calculation][planning][sqe]") {
    for (const auto method_id : {"sqe", "sqeq0"}) {
        auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules =
                core::MoleculeCollection{std::vector{make_fractional_sqe_source_molecule()}},
            .parameter_sets = {make_sqe_embedding_parameters(method_id)},
            .method_id = method_id,
            .fixed_charge_embedding = calculation::FixedChargeEmbedding{
                .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 0.25},
                            {.molecule_index = 0, .atom_index = 3, .charge = -0.1}}}});
        REQUIRE(assessment.plans().size() == 1);
        const auto result = calculation::calculate(assessment);
        REQUIRE(result.calculated());
        const auto& charges = result.charges->assignment(0).charges;
        CHECK(charges[1] == 0.25);
        CHECK(charges[3] == -0.1);
        CHECK(std::abs(charges[0] + charges[2]) < 1.0e-12);
        CHECK(std::abs(charges[0] + charges[1] + charges[2] + charges[3] - 0.15) < 1.0e-12);
    }
}

TEST_CASE("unsupported methods reject embeddings before parameter classification",
          "[calculation][planning]") {
    for (const auto method_id : {"formal", "peoe"}) {
        auto request = calculation::AssessmentRequest{
            .molecules =
                core::MoleculeCollection{std::vector{make_interleaved_hydrogen_magnesium_oxygen()}},
            .method_id = method_id,
            .fixed_charge_embedding = calculation::FixedChargeEmbedding{
                .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}}}};
        const auto assessment = calculation::assess(std::move(request));

        CHECK(assessment.plans().empty());
        REQUIRE(assessment.rejections().size() == 1);
        CHECK(assessment.rejections()[0].method_id == method_id);
        CHECK_FALSE(assessment.rejections()[0].policy.has_value());
        REQUIRE(assessment.rejections()[0].issues.size() == 1);
        CHECK(std::get<methods::PrerequisiteIssue>(assessment.rejections()[0].issues[0]).kind ==
              methods::PrerequisiteIssueKind::unsupported_embedding);
    }
}

TEST_CASE("embedded EEM permits full execution and blocks reduced modes",
          "[calculation][planning]") {
    const auto assess_mode = [](const calculation::ExecutionSelectionKind kind) {
        return calculation::assess(calculation::AssessmentRequest{
            .molecules =
                core::MoleculeCollection{std::vector{make_interleaved_hydrogen_magnesium_oxygen()}},
            .parameter_sets = {make_eem_parameters()},
            .method_id = "eem",
            .execution_selection =
                calculation::ExecutionSelection{
                    kind, kind == calculation::ExecutionSelectionKind::full
                              ? std::optional<double>{}
                              : std::optional{calculation::minimum_reduced_radius}},
            .resource_policy = {.cutoff_atom_threshold = 0, .cover_atom_threshold = 0},
            .fixed_charge_embedding = calculation::FixedChargeEmbedding{
                .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}}}});
    };
    for (const auto& [selection, mode] : {std::pair{calculation::ExecutionSelectionKind::cutoff,
                                                    calculation::ExecutionMode::cutoff},
                                          std::pair{calculation::ExecutionSelectionKind::cover,
                                                    calculation::ExecutionMode::cover}}) {
        const auto assessment = assess_mode(selection);
        CHECK(assessment.plans().empty());
        REQUIRE(assessment.rejections().size() == 1);
        REQUIRE(assessment.rejections()[0].policy.has_value());
        CHECK(assessment.rejections()[0].policy->mode() == mode);
        REQUIRE(assessment.rejections()[0].issues.size() == 1);
        const auto& issue = std::get<methods::ExecutionIssue>(assessment.rejections()[0].issues[0]);
        CHECK(issue.kind == methods::ExecutionIssueKind::unsupported_execution_mode);
        CHECK(issue.message.contains("fixed charge embedding supports full execution only"));
    }

    const auto full_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules =
            core::MoleculeCollection{std::vector{make_interleaved_hydrogen_magnesium_oxygen()}},
        .parameter_sets = {make_eem_parameters()},
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}}}});
    REQUIRE(full_assessment.plans().size() == 1);
    CHECK(full_assessment.default_plan()->policy().mode() == calculation::ExecutionMode::full);
    CHECK(calculation::calculate(full_assessment).calculated());

    auto automatic_request = calculation::AssessmentRequest{
        .molecules =
            core::MoleculeCollection{std::vector{make_interleaved_hydrogen_magnesium_oxygen()}},
        .parameter_sets = {make_eem_parameters()},
        .resource_policy = {.cutoff_atom_threshold = 0, .cover_atom_threshold = 0},
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}}}};
    const auto automatic = calculation::assess(std::move(automatic_request));
    CHECK(automatic.plans().empty());
    auto eem_execution_rejections = 0;
    for (const auto& rejection : automatic.rejections()) {
        if (rejection.method_id == "eem") {
            ++eem_execution_rejections;
            REQUIRE(rejection.policy.has_value());
            REQUIRE(rejection.issues.size() == 1);
            const auto& issue = std::get<methods::ExecutionIssue>(rejection.issues[0]);
            CHECK(issue.kind == (rejection.policy->mode() == calculation::ExecutionMode::full
                                     ? methods::ExecutionIssueKind::resource_threshold_exceeded
                                     : methods::ExecutionIssueKind::unsupported_execution_mode));
        } else {
            CHECK_FALSE(rejection.policy.has_value());
            REQUIRE(rejection.issues.size() == 1);
            const auto& issue = std::get<methods::PrerequisiteIssue>(rejection.issues[0]);
            CHECK(issue.kind == ((rejection.method_id == "sqe" || rejection.method_id == "sqeq0")
                                     ? methods::PrerequisiteIssueKind::missing_parameters
                                     : methods::PrerequisiteIssueKind::unsupported_embedding));
        }
    }
    CHECK(eem_execution_rejections == 3);
}

TEST_CASE("embedded assessment retains original molecules and rejections after moves",
          "[calculation][planning]") {
    const auto create_assessment = [] {
        auto request = calculation::AssessmentRequest{
            .molecules =
                core::MoleculeCollection{std::vector{make_interleaved_hydrogen_magnesium_oxygen()}},
            .parameter_sets = {make_eem_parameters()},
            .method_id = "eem",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
            .fixed_charge_embedding = calculation::FixedChargeEmbedding{
                .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}},
                .charge_provenance = "owned sources"}};
        return calculation::assess(std::move(request));
    };
    auto assessment = create_assessment();
    auto moved_assessment = std::move(assessment);

    REQUIRE(moved_assessment.molecules().size() == 1);
    CHECK(moved_assessment.molecules()[0].atom_count() == 3);
    CHECK(moved_assessment.molecules()[0].atom(1).name() == "fixed-Mg");
    REQUIRE(moved_assessment.plans().size() == 1);
    CHECK(moved_assessment.plans()[0].policy().mode() == calculation::ExecutionMode::full);
    const auto moved_result = calculation::calculate(moved_assessment);
    REQUIRE(moved_result.calculated());
    REQUIRE(moved_result.effective.has_value());
    REQUIRE(moved_result.effective->fixed_charge_embedding.has_value());
}

TEST_CASE("embedded parameter rejections use original atom indices and descriptions",
          "[calculation][planning]") {
    const auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules =
            core::MoleculeCollection{std::vector{make_interleaved_hydrogen_magnesium_oxygen()}},
        .parameter_sets = {make_hydrogen_only_eem_parameters()},
        .method_id = "eem",
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}}}});

    CHECK(assessment.plans().empty());
    REQUIRE(assessment.rejections().size() == 1);
    CHECK(assessment.rejections()[0].parameter_set_id == "test-eem-h-only");
    REQUIRE(assessment.rejections()[0].issues.size() == 1);
    const auto& issue = std::get<methods::PrerequisiteIssue>(assessment.rejections()[0].issues[0]);
    CHECK(issue.kind == methods::PrerequisiteIssueKind::parameter_classification_failed);
    CHECK(issue.molecule_index == 0);
    CHECK(issue.atom_index == 2);
    CHECK(issue.message.contains("parameter set 'test-eem-h-only'"));
    CHECK(issue.message.contains("atom 3 (source name 'active-O', O, formal charge 0)"));
}

TEST_CASE("embedded geometry diagnostics label active-subsystem atom numbering",
          "[calculation][planning]") {
    const auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{
            make_interleaved_hydrogen_magnesium_oxygen(core::Position{0.0, 0.0, 0.0})}},
        .parameter_sets = {make_eem_parameters()},
        .method_id = "eem",
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}}}});

    CHECK(assessment.plans().empty());
    REQUIRE(assessment.rejections().size() == 1);
    REQUIRE(assessment.rejections()[0].issues.size() == 1);
    const auto& issue = std::get<methods::PrerequisiteIssue>(assessment.rejections()[0].issues[0]);
    CHECK(issue.kind == methods::PrerequisiteIssueKind::invalid_geometry);
    CHECK(issue.atom_index == 2);
    CHECK(issue.message.contains("molecule 1 ('interleaved-H-Mg-O')"));
    CHECK(issue.message.contains("active-subsystem atom numbering"));
}

TEST_CASE("fixed charge embedding source selectors are validated", "[calculation][planning]") {
    const auto check_invalid = [](std::vector<core::Molecule> molecules,
                                  std::vector<calculation::FixedAtomCharge> sources,
                                  const std::string_view diagnostic) {
        auto request = calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::move(molecules)},
            .fixed_charge_embedding =
                calculation::FixedChargeEmbedding{.sources = std::move(sources)}};
        auto rejected = false;
        try {
            static_cast<void>(calculation::assess(std::move(request)));
        } catch (const std::invalid_argument& error) {
            rejected = true;
            CHECK(std::string_view{error.what()}.contains(diagnostic));
        }
        CHECK(rejected);
    };

    const auto ion_pair = make_isolated_ion_pair();
    check_invalid({ion_pair}, {{2, 0, 1.0}}, "molecule index 2 is out of range");
    check_invalid({ion_pair}, {{std::numeric_limits<std::size_t>::max(), 0, 1.0}},
                  "is out of range (molecule count 1)");
    check_invalid({ion_pair}, {{0, std::numeric_limits<std::size_t>::max(), 1.0}},
                  "out of range (atom count 2)");
    check_invalid({ion_pair}, {{0, 0, 1.0}, {0, 0, -1.0}}, "duplicate fixed charge source");
    for (const auto charge :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity()}) {
        check_invalid({ion_pair}, {{0, 0, charge}}, "has non-finite charge");
    }
    check_invalid({make_bonded_pair()}, {{0, 0, 0.5}}, "is bonded to atom 1");
    check_invalid({make_bonded_pair()}, {{0, 1, 0.5}}, "is bonded to atom 0");
    check_invalid({make_bonded_pair()}, {{0, 0, 0.5}, {0, 1, -0.5}}, "is bonded to atom 1");
    check_invalid({core::Molecule{std::vector{core::Atom{6}}}, ion_pair}, {{0, 0, 0.0}},
                  "leaves no active atoms in molecule 0");
}

TEST_CASE("fixed charge embedding validates geometry only for affected molecules",
          "[calculation][planning]") {
    const auto make_request = [](core::MoleculeCollection molecules,
                                 std::vector<calculation::FixedAtomCharge> sources) {
        return calculation::AssessmentRequest{
            .molecules = std::move(molecules),
            .fixed_charge_embedding =
                calculation::FixedChargeEmbedding{.sources = std::move(sources)}};
    };
    const auto pair_atoms = std::vector{core::Atom{6}, core::Atom{6}};
    const auto separated_pair =
        std::vector<std::vector<core::Position>>{{{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}}};

    auto missing_conformer = make_request(
        core::MoleculeCollection{std::vector{make_molecule_with_positions(pair_atoms, {})}},
        {{0, 0, 1.0}});
    CHECK(std::string_view{assessment_error(std::move(missing_conformer))}.contains(
        "requires a conformer in molecule 0"));

    for (std::size_t atom_index = 0; atom_index < 2; ++atom_index) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            auto later_positions = separated_pair.front();
            auto& position = later_positions[atom_index];
            const auto value = axis == 0   ? std::numeric_limits<double>::quiet_NaN()
                               : axis == 1 ? std::numeric_limits<double>::infinity()
                                           : -std::numeric_limits<double>::infinity();
            if (axis == 0) {
                position.x = value;
            } else if (axis == 1) {
                position.y = value;
            } else {
                position.z = value;
            }
            auto nonfinite = make_request(
                core::MoleculeCollection{std::vector{make_molecule_with_positions(
                    pair_atoms, {separated_pair.front(), std::move(later_positions)})}},
                {{0, 0, 1.0}});
            CHECK(std::string_view{assessment_error(std::move(nonfinite))}.contains(
                "non-finite coordinates at molecule 0, conformer 1, atom " +
                std::to_string(atom_index)));
        }
    }

    auto coincident_later = make_request(
        core::MoleculeCollection{std::vector{make_molecule_with_positions(
            pair_atoms, {separated_pair.front(), {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}})}},
        {{0, 0, 0.0}});
    CHECK(std::string_view{assessment_error(std::move(coincident_later))}.contains(
        "molecule 0, conformer 1, atom 0 coincides with active atom 1"));

    auto swapped_conformers = make_request(
        core::MoleculeCollection{std::vector{make_molecule_with_positions(
            pair_atoms, {separated_pair.front(), {{2.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}})}},
        {{0, 0, 0.0}});
    const auto swapped_assessment = calculation::assess(std::move(swapped_conformers));
    CHECK(swapped_assessment.plans().empty());
    CHECK_FALSE(swapped_assessment.rejections().empty());

    auto near_nonzero =
        make_request(core::MoleculeCollection{std::vector{make_molecule_with_positions(
                         pair_atoms, {{{0.0, 0.0, 0.0}, {1e-200, 0.0, 0.0}}})}},
                     {{0, 0, 0.0}});
    const auto near_assessment = calculation::assess(std::move(near_nonzero));
    CHECK(near_assessment.plans().empty());
    CHECK_FALSE(near_assessment.rejections().empty());

    auto coincident_sources =
        make_request(core::MoleculeCollection{std::vector{make_molecule_with_positions(
                         std::vector{core::Atom{6}, core::Atom{6}, core::Atom{6}},
                         {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}}})}},
                     {{0, 0, 0.5}, {0, 1, -0.25}});
    const auto sources_assessment = calculation::assess(std::move(coincident_sources));
    CHECK(sources_assessment.plans().empty());
    CHECK_FALSE(sources_assessment.rejections().empty());

    auto unrelated = make_request(
        core::MoleculeCollection{std::vector{
            make_molecule_with_positions(pair_atoms, separated_pair),
            make_molecule_with_positions(pair_atoms, {{{0.0, 0.0, 0.0}, {3.0, 0.0, 0.0}}}),
            make_molecule_with_positions(pair_atoms, {}),
            make_molecule_with_positions(
                pair_atoms,
                {{{std::numeric_limits<double>::infinity(), 0.0, 0.0}, {4.0, 0.0, 0.0}}})}},
        {{0, 0, 1.0}});
    const auto unrelated_assessment = calculation::assess(std::move(unrelated));
    CHECK(unrelated_assessment.plans().empty());
    CHECK_FALSE(unrelated_assessment.rejections().empty());
}

TEST_CASE("fixed charge embedding records active and original totals independently",
          "[calculation][planning]") {
    auto audit_vs_model_request = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_isolated_ion_pair()}},
        .parameter_sets = {make_eem_parameters()},
        .method_id = "eem",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = -0.5}}}};
    const auto assessment = calculation::assess(audit_vs_model_request);
    REQUIRE(assessment.plans().size() == 1);
    REQUIRE(audit_vs_model_request.fixed_charge_embedding.has_value());
    CHECK(audit_vs_model_request.fixed_charge_embedding->sources[0].charge == -0.5);
    CHECK(core::total_formal_charge(audit_vs_model_request.molecules[0]) == 2.0);
    const auto result = calculation::calculate(assessment);
    REQUIRE(result.calculated());
    REQUIRE(result.effective->fixed_charge_embedding.has_value());
    CHECK(result.effective->fixed_charge_embedding->charge_totals[0].original_total_charge == 2.0);
    CHECK(result.effective->fixed_charge_embedding->charge_totals[0].active_total_charge == 0.0);
    CHECK(result.charges->assignment(0).charges[1] == -0.5);
    CHECK(std::abs(result.charges->assignment(0).charges.total() + 0.5) < 1.0e-12);
}

TEST_CASE("empty fixed charge embedding is equivalent to absence", "[calculation][planning]") {
    auto make_request = [](const bool with_empty_embedding) {
        auto request = calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {make_eem_parameters()}};
        if (with_empty_embedding) {
            request.fixed_charge_embedding =
                calculation::FixedChargeEmbedding{.charge_provenance = "ignored when empty"};
        }
        return request;
    };

    const auto ordinary = calculation::assess(make_request(false));
    const auto empty = calculation::assess(make_request(true));
    REQUIRE_FALSE(ordinary.plans().empty());
    REQUIRE(empty.plans().size() == ordinary.plans().size());
    for (std::size_t index = 0; index < ordinary.plans().size(); ++index) {
        CHECK(empty.plans()[index].candidate().method->id() ==
              ordinary.plans()[index].candidate().method->id());
        CHECK(empty.plans()[index].policy().mode() == ordinary.plans()[index].policy().mode());
        CHECK(empty.plans()[index].policy().radius() == ordinary.plans()[index].policy().radius());
    }
    const auto ordinary_result = calculation::calculate(ordinary);
    const auto empty_result = calculation::calculate(empty);
    REQUIRE(ordinary_result.calculated());
    REQUIRE(empty_result.calculated());
    CHECK_FALSE(ordinary_result.effective->fixed_charge_embedding.has_value());
    CHECK_FALSE(empty_result.effective->fixed_charge_embedding.has_value());
    CHECK(std::ranges::equal(empty_result.charges->assignment(0).charges.values(),
                             ordinary_result.charges->assignment(0).charges.values()));
}
