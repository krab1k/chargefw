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

    CHECK(assessment.plans().empty());
    REQUIRE(assessment.molecules().size() == 3);
    CHECK(assessment.molecules()[1].atom(1).atomic_number() == 12);
    CHECK(assessment.molecules()[1].atom(1).name() == "fixed-Mg");
    REQUIRE(assessment.rejections().size() == 3);
    for (const auto& rejection : assessment.rejections()) {
        CHECK(rejection.method_id == "eem");
        REQUIRE(rejection.policy.has_value());
        REQUIRE(rejection.issues.size() == 1);
        CHECK(std::get<methods::ExecutionIssue>(rejection.issues[0]).kind ==
              methods::ExecutionIssueKind::unsupported_execution_mode);
    }
    const auto execution = calculation::calculate(assessment);
    CHECK(execution.status == calculation::ExecutionStatus::no_executable_plan);
    CHECK_FALSE(execution.charges.has_value());
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

TEST_CASE("embedded EEM assessments block every execution mode without resource overrides",
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
    for (const auto& [selection, mode] :
         {std::pair{calculation::ExecutionSelectionKind::full, calculation::ExecutionMode::full},
          std::pair{calculation::ExecutionSelectionKind::cutoff,
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
        CHECK(issue.message.contains("fixed charge embedding execution is not connected yet"));
    }

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
            CHECK(std::get<methods::ExecutionIssue>(rejection.issues[0]).kind ==
                  methods::ExecutionIssueKind::unsupported_execution_mode);
        } else {
            CHECK_FALSE(rejection.policy.has_value());
            REQUIRE(rejection.issues.size() == 1);
            CHECK(std::get<methods::PrerequisiteIssue>(rejection.issues[0]).kind ==
                  methods::PrerequisiteIssueKind::unsupported_embedding);
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
    CHECK(moved_assessment.plans().empty());
    REQUIRE(moved_assessment.rejections().size() == 1);
    CHECK(moved_assessment.rejections()[0].method_id == "eem");
    REQUIRE(moved_assessment.rejections()[0].policy.has_value());
    CHECK(moved_assessment.rejections()[0].policy->mode() == calculation::ExecutionMode::full);
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

TEST_CASE("fixed charge embedding validates per-molecule charge budgets",
          "[calculation][planning]") {
    const auto make_request = [](core::MoleculeCollection molecules,
                                 std::vector<calculation::FixedAtomCharge> sources) {
        return calculation::AssessmentRequest{
            .molecules = std::move(molecules),
            .fixed_charge_embedding =
                calculation::FixedChargeEmbedding{.sources = std::move(sources)}};
    };
    const auto three_atoms = std::vector{core::Atom{6}, core::Atom{6}, core::Atom{6}};
    const auto pair_atoms = std::vector{core::Atom{6}, core::Atom{6}};
    const auto three_positions = std::vector<std::vector<core::Position>>{
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}}};
    const auto maximum = std::numeric_limits<double>::max();
    auto overflowing_sum =
        make_request(core::MoleculeCollection{std::vector{
                         make_molecule_with_positions(three_atoms, three_positions)}},
                     {{0, 0, maximum}, {0, 1, maximum}});
    CHECK(std::string_view{assessment_error(std::move(overflowing_sum))}.contains(
        "fixed charge source sum is non-finite in molecule 0"));

    auto separate_targets = make_request(
        core::MoleculeCollection{std::vector{
            make_molecule_with_positions(pair_atoms, {{{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}}}),
            make_molecule_with_positions(pair_atoms, {{{0.0, 0.0, 0.0}, {3.0, 0.0, 0.0}}})}},
        {{0, 0, maximum}, {1, 0, maximum}});
    const auto separate_assessment = calculation::assess(std::move(separate_targets));
    CHECK(separate_assessment.plans().empty());
    CHECK_FALSE(separate_assessment.rejections().empty());

    auto mismatched_budget = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_isolated_ion_pair()}},
        .parameter_sets = {make_eem_parameters()},
        .method_id = "eem",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
        .fixed_charge_embedding = calculation::FixedChargeEmbedding{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = -0.5}}}};
    const auto mismatched_assessment = calculation::assess(mismatched_budget);
    CHECK(mismatched_assessment.plans().empty());
    REQUIRE(mismatched_assessment.rejections().size() == 1);
    REQUIRE(mismatched_assessment.rejections()[0].policy.has_value());
    CHECK(mismatched_assessment.rejections()[0].policy->mode() == calculation::ExecutionMode::full);
    REQUIRE(mismatched_budget.fixed_charge_embedding.has_value());
    CHECK(mismatched_budget.fixed_charge_embedding->sources[0].charge == -0.5);
    CHECK(core::total_formal_charge(mismatched_budget.molecules[0]) == 2.0);
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
    CHECK(std::ranges::equal(empty_result.charges->assignment(0).charges.values(),
                             ordinary_result.charges->assignment(0).charges.values()));
}
