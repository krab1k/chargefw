#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cmath>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

[[nodiscard]] auto calculate_application(calculation::AssessmentRequest request)
    -> calculation::ExecutionResult {
    const auto max_threads = request.resource_policy.max_threads;
    auto assessment = calculation::assess(std::move(request));
    return calculation::calculate(std::move(assessment), max_threads);
}

class ZeroFragmentMethod final : public methods::Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{.id = "zero-fragment",
                                                          .name = "Zero fragment",
                                                          .full_name = "Zero fragment",
                                                          .publication = std::nullopt,
                                                          .priority = 0};
        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        auto requirements = methods::MethodRequirements{};
        requirements.coordinates = true;
        requirements.resources.supports_cutoff = true;
        requirements.resources.supports_cover = true;
        requirements.resources.fragment_target_charge_policy =
            methods::FragmentTargetChargePolicy::proportional_to_atom_count;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> charges::AtomicCharges override {
        return charges::AtomicCharges{std::vector<double>(input.molecule().atom_count(), 0.0)};
    }
};

auto make_eem_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "test-eem", .method_id = "eem", .name = "Test EEM"},
        parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
        parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 10.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 10.0}}}}}};
}

auto make_qeq_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "test-qeq", .method_id = "qeq", .name = "Test QEq"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 4.5280},
                                                    {.name = "hardness", .value = 13.8904}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 8.741},
                                                    {.name = "hardness", .value = 13.364}}}}}};
}

auto make_invalid_qeq_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "invalid-qeq", .method_id = "qeq", .name = "Invalid QEq"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 4.5280},
                                                    {.name = "hardness", .value = 0.0}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 8.741},
                                                    {.name = "hardness", .value = 13.364}}}}}};
}

auto make_eqeqc_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-eqeqc", .method_id = "eqeqc", .name = "Test EQeq+C"},
        parameters::CommonParameters{{{.name = "alpha", .value = 1.0}}},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "Dz", .value = 0.1}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "Dz", .value = 0.2}}}}}};
}

auto make_abeem_parameters() -> parameters::ParameterSet {
    const auto bond_key = parameters::BondParameterKey{
        .first_atom = chargefw::test::plain_atom_key(8),
        .second_atom = chargefw::test::plain_atom_key(1),
        .bond = {.classification = parameters::BondParameterClassificationKind::PLAIN,
                 .type = "*"}};
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-abeem", .method_id = "abeem", .name = "Test ABEEM"},
        parameters::CommonParameters{{{.name = "k", .value = 1.0}}},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "a", .value = 1.0},
                                                    {.name = "b", .value = 10.0},
                                                    {.name = "c", .value = 0.5}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "a", .value = 2.0},
                                                    {.name = "b", .value = 10.0},
                                                    {.name = "c", .value = 0.5}}}}},
        parameters::BondParameters{{{.key = bond_key,
                                     .parameters = {{.name = "A", .value = 1.0},
                                                    {.name = "B", .value = 10.0},
                                                    {.name = "C", .value = 0.5},
                                                    {.name = "D", .value = 0.5}}}}}};
}

auto make_sqe_parameters(const std::string_view method_id, const bool parameterized_initial_charge,
                         const bool zero_widths = false) -> parameters::ParameterSet {
    auto hydrogen_parameters = std::vector<parameters::NamedParameter>{
        {.name = "electronegativity", .value = 4.5280},
        {.name = "hardness", .value = 13.8904},
        {.name = "width", .value = zero_widths ? 0.0 : 1.0}};
    auto oxygen_parameters = std::vector<parameters::NamedParameter>{
        {.name = "electronegativity", .value = 8.741},
        {.name = "hardness", .value = 13.364},
        {.name = "width", .value = zero_widths ? 0.0 : 1.0}};
    if (parameterized_initial_charge) {
        hydrogen_parameters.push_back({.name = "q0", .value = 0.25});
        oxygen_parameters.push_back({.name = "q0", .value = -0.5});
    }

    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "test-" + std::string{method_id},
                                         .method_id = std::string{method_id},
                                         .name = "Test SQE-family parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = std::move(hydrogen_parameters)},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = std::move(oxygen_parameters)}}},
        parameters::BondParameters{{{.key = chargefw::test::single_bond_key(1, 8),
                                     .parameters = {{.name = "kappa", .value = 1.0}}}}}};
}

auto make_charged_water() -> core::Molecule {
    return core::Molecule{
        std::vector{core::Atom{8, 1, "O"}, core::Atom{1, 0, "H1"}, core::Atom{1, 0, "H2"}},
        std::vector{core::Bond{0, 1, core::BondOrder::SINGLE},
                    core::Bond{0, 2, core::BondOrder::SINGLE}},
        std::vector{core::Conformer{{core::Position{.x = 0.0, .y = 0.0, .z = 0.0},
                                     core::Position{.x = 0.9572, .y = 0.0, .z = 0.0},
                                     core::Position{.x = -0.2390, .y = 0.9270, .z = 0.0}}}},
        "charged-water"};
}

auto assert_reduced_matches_full(
    const std::string_view method_id,
    const std::vector<parameters::ParameterSet>& parameter_sets = {},
    core::Molecule molecule = chargefw::test::make_two_conformer_water()) -> void {
    const auto molecules = core::MoleculeCollection{std::vector{molecule, std::move(molecule)}};
    const auto full = calculate_application(
        calculation::AssessmentRequest{.molecules = molecules,
                                       .parameter_sets = parameter_sets,
                                       .method_id = std::string{method_id},
                                       .resource_policy = {.max_threads = 2}});
    CHECK(full.calculated());

    for (const auto selection_kind : {calculation::ExecutionSelectionKind::cutoff,
                                      calculation::ExecutionSelectionKind::cover}) {
        const auto reduced = calculate_application(calculation::AssessmentRequest{
            .molecules = molecules,
            .parameter_sets = parameter_sets,
            .method_id = std::string{method_id},
            .execution_selection = calculation::ExecutionSelection{selection_kind, 8.0},
            .resource_policy = {.max_threads = 2}});

        CHECK(reduced.calculated());
        CHECK(reduced.execution_policy.has_value());
        CHECK(reduced.execution_policy->mode() ==
              (selection_kind == calculation::ExecutionSelectionKind::cutoff
                   ? calculation::ExecutionMode::cutoff
                   : calculation::ExecutionMode::cover));
        CHECK(reduced.execution_policy->charge_correction() ==
              calculation::ChargeCorrectionPolicy::uniform);
        CHECK(full.charges->size() == reduced.charges->size());

        for (std::size_t assignment_index = 0; assignment_index < full.charges->size();
             ++assignment_index) {
            const auto& full_charges = full.charges->assignment(assignment_index).charges;
            const auto& reduced_charges = reduced.charges->assignment(assignment_index).charges;
            CHECK(full_charges.size() == reduced_charges.size());
            for (std::size_t atom_index = 0; atom_index < full_charges.size(); ++atom_index) {
                CHECK(std::abs(full_charges[atom_index] - reduced_charges[atom_index]) < 1.0e-10);
            }
        }
    }
}

} // namespace

TEST_CASE("reduced execution validates inputs, correction, and mode selection",
          "[calculation][reduced-execution]") {
    const ZeroFragmentMethod zero_method;
    const auto charged_molecule = core::Molecule{
        std::vector{core::Atom{1, 1}, core::Atom{1, 0}},
        {},
        std::vector{core::Conformer{{core::Position{.x = 0.0}, core::Position{.x = 20.0}}}},
        "charged"};
    const auto collection = core::MoleculeCollection{std::vector{charged_molecule}};
    const features::PreparedMoleculeCollection prepared{collection};
    const methods::ApplicableMethod selected{.method = &zero_method, .parameter_set = nullptr};

    const methods::ApplicableMethod invalid_selected{
        .method = &zero_method, .parameter_set = nullptr, .classifications = {{}}};
    for (const auto mode : {calculation::ExecutionMode::full, calculation::ExecutionMode::cutoff,
                            calculation::ExecutionMode::cover}) {
        const auto policy = mode == calculation::ExecutionMode::full
                                ? calculation::ExecutionPolicy{}
                                : calculation::ExecutionPolicy{mode, 8.0};
        const auto calculate_with_invalid_classification = [&] -> void {
            static_cast<void>(calculation::calculate(
                {.molecules = prepared, .selected = invalid_selected, .execution_policy = policy}));
        };
        CHECK_THROWS_AS(calculate_with_invalid_classification(), std::invalid_argument);
    }

    const auto no_conformer_collection = core::MoleculeCollection{
        std::vector{core::Molecule{std::vector{core::Atom{1}}, {}, {}, "no-conformer"}}};
    const features::PreparedMoleculeCollection no_conformer_prepared{no_conformer_collection};
    for (const auto mode : {calculation::ExecutionMode::full, calculation::ExecutionMode::cutoff,
                            calculation::ExecutionMode::cover}) {
        const auto policy = mode == calculation::ExecutionMode::full
                                ? calculation::ExecutionPolicy{}
                                : calculation::ExecutionPolicy{mode, 8.0};
        const calculation::CalculationRequest request{
            .molecules = no_conformer_prepared, .selected = selected, .execution_policy = policy};
        const auto calculate_without_conformer = [&request] -> void {
            static_cast<void>(calculation::calculate(request));
        };
        CHECK_THROWS_AS(calculate_without_conformer(), std::invalid_argument);
    }

    const auto corrected = calculation::calculate(
        {.molecules = prepared,
         .selected = selected,
         .execution_policy =
             calculation::ExecutionPolicy{calculation::ExecutionMode::cutoff, 8.0,
                                          calculation::ChargeCorrectionPolicy::uniform},
         .max_threads = 2});
    CHECK(corrected.charges.assignment(0).charges[0] == 0.5);
    CHECK(corrected.charges.assignment(0).charges[1] == 0.5);

    const auto uncorrected = calculation::calculate(
        {.molecules = prepared,
         .selected = selected,
         .execution_policy =
             calculation::ExecutionPolicy{calculation::ExecutionMode::cutoff, 8.0,
                                          calculation::ChargeCorrectionPolicy::none},
         .max_threads = 2});
    CHECK(uncorrected.charges.assignment(0).charges[0] == 0.0);
    CHECK(uncorrected.charges.assignment(0).charges[1] == 0.0);

    const auto cover_corrected =
        calculation::calculate({.molecules = prepared,
                                .selected = selected,
                                .execution_policy = calculation::ExecutionPolicy{
                                    calculation::ExecutionMode::cover, 8.0,
                                    calculation::ChargeCorrectionPolicy::uniform}});
    CHECK(cover_corrected.charges.assignment(0).charges[0] == 0.5);
    CHECK(cover_corrected.charges.assignment(0).charges[1] == 0.5);

    assert_reduced_matches_full("eem", {make_eem_parameters()});

    const auto automatic_cutoff = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {make_eem_parameters()},
        .method_id = "eem",
        .resource_policy = {.cutoff_atom_threshold = 2}});
    CHECK(automatic_cutoff.calculated());
    CHECK(automatic_cutoff.execution_policy->mode() == calculation::ExecutionMode::cutoff);
    CHECK(automatic_cutoff.execution_policy->radius() ==
          std::optional<double>{calculation::default_automatic_reduced_radius});

    const auto overridden_automatic_cutoff = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {make_eem_parameters()},
        .method_id = "eem",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::automatic, 8.0},
        .resource_policy = {.cutoff_atom_threshold = 2}});
    CHECK(overridden_automatic_cutoff.calculated());
    CHECK(overridden_automatic_cutoff.execution_policy->mode() ==
          calculation::ExecutionMode::cutoff);
    CHECK(overridden_automatic_cutoff.execution_policy->radius() == std::optional<double>{8.0});

    const auto automatic_cover = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {make_eem_parameters()},
        .method_id = "eem",
        .resource_policy = {.cutoff_atom_threshold = 2, .cover_atom_threshold = 2}});
    CHECK(automatic_cover.calculated());
    CHECK(automatic_cover.execution_policy->mode() == calculation::ExecutionMode::cover);

    const auto explicit_cutoff_above_cover_threshold =
        calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {make_eem_parameters()},
            .method_id = "eem",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff,
                                                calculation::minimum_reduced_radius},
            .resource_policy = {.cutoff_atom_threshold = 2, .cover_atom_threshold = 2}});
    CHECK(explicit_cutoff_above_cover_threshold.calculated());
    CHECK(explicit_cutoff_above_cover_threshold.execution_policy->mode() ==
          calculation::ExecutionMode::cutoff);
    CHECK(explicit_cutoff_above_cover_threshold.execution_issues.size() == 1);

    assert_reduced_matches_full("qeq", {make_qeq_parameters()});
    assert_reduced_matches_full("eqeq");
    assert_reduced_matches_full("eqeqc", {make_eqeqc_parameters()});
    assert_reduced_matches_full("abeem", {make_abeem_parameters()});
    assert_reduced_matches_full("sqe", {make_sqe_parameters("sqe", false)});
    assert_reduced_matches_full("sqeq0", {make_sqe_parameters("sqeq0", false)});
    assert_reduced_matches_full("sqeqp", {make_sqe_parameters("sqeqp", true)});
    assert_reduced_matches_full("sqe", {make_sqe_parameters("sqe", false, true)});
    assert_reduced_matches_full("sqeq0", {make_sqe_parameters("sqeq0", false, true)});
    assert_reduced_matches_full("sqeqp", {make_sqe_parameters("sqeqp", true, true)});
    assert_reduced_matches_full("sqe", {make_sqe_parameters("sqe", false)}, make_charged_water());
    assert_reduced_matches_full("sqeq0", {make_sqe_parameters("sqeq0", false)},
                                make_charged_water());
    assert_reduced_matches_full("sqeqp", {make_sqe_parameters("sqeqp", true)},
                                make_charged_water());
}

TEST_CASE("reduced solver failures retain method and target context",
          "[calculation][reduced-execution]") {
    for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                            calculation::ExecutionSelectionKind::cover}) {
        const auto calculate_invalid_qeq = [mode] {
            static_cast<void>(calculate_application(calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
                .parameter_sets = {make_invalid_qeq_parameters()},
                .method_id = "qeq",
                .parameter_set_id = "invalid-qeq",
                .execution_selection =
                    calculation::ExecutionSelection{mode, calculation::minimum_reduced_radius},
                .resource_policy = {.max_threads = 1}}));
        };

        try {
            calculate_invalid_qeq();
            CHECK(false);
        } catch (const std::runtime_error& error) {
            const auto message = std::string_view{error.what()};
            CHECK(message.contains(mode == calculation::ExecutionSelectionKind::cutoff
                                       ? "center atom"
                                       : "pivot atom"));
            CHECK(message.contains("method 'qeq'"));
            CHECK(message.contains("molecule 'water'"));
            CHECK(message.contains("conformer 0"));
        }
    }
}

TEST_CASE("reduced execution preserves mixed source target order",
          "[calculation][reduced-execution]") {
    const auto collection = core::MoleculeCollection{
        std::vector{chargefw::test::make_two_conformer_water(), chargefw::test::make_water()},
        "mixed-water"};
    const features::PreparedMoleculeCollection prepared{collection};
    const ZeroFragmentMethod method;
    const methods::ApplicableMethod selected{.method = &method, .parameter_set = nullptr};

    for (const auto mode :
         {calculation::ExecutionMode::cutoff, calculation::ExecutionMode::cover}) {
        const auto result = calculation::calculate(
            {.molecules = prepared,
             .selected = selected,
             .execution_policy =
                 calculation::ExecutionPolicy{mode, calculation::minimum_reduced_radius,
                                              calculation::ChargeCorrectionPolicy::uniform},
             .max_threads = 2});

        REQUIRE(result.charges.size() == 3);
        CHECK(result.charges.assignment(0).target.molecule_index == 0);
        CHECK(result.charges.assignment(0).target.conformer_index == std::optional<std::size_t>{0});
        CHECK(result.charges.assignment(1).target.molecule_index == 0);
        CHECK(result.charges.assignment(1).target.conformer_index == std::optional<std::size_t>{1});
        CHECK(result.charges.assignment(2).target.molecule_index == 1);
        CHECK(result.charges.assignment(2).target.conformer_index == std::optional<std::size_t>{0});
    }
}
