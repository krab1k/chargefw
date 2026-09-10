#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/features/spatial_fragment.h>
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

#include <algorithm>
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
    auto assessment = calculation::assess(std::move(request));
    return calculation::calculate(assessment, 1);
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
        requirements.resources.reduced_charge_policy =
            methods::ReducedChargePolicy::uniform_target_global;
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

auto make_sfkeem_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-sfkeem", .method_id = "sfkeem", .name = "Test SFKEEM"},
        parameters::CommonParameters{{{.name = "sigma", .value = 1.0}}},
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

auto make_extended_components(const int isolated_formal_charge) -> core::Molecule {
    constexpr auto chain_atom_count = std::size_t{12};
    auto atoms = std::vector<core::Atom>{};
    auto bonds = std::vector<core::Bond>{};
    auto positions = std::vector<core::Position>{};
    atoms.reserve(chain_atom_count + 1);
    bonds.reserve(chain_atom_count - 1);
    positions.reserve(chain_atom_count + 1);

    for (std::size_t atom_index = 0; atom_index < chain_atom_count; ++atom_index) {
        const auto atomic_number = atom_index % 2 == 0 ? 8 : 1;
        atoms.emplace_back(atomic_number, 0, "chain");
        positions.push_back(core::Position{.x = 1.4 * static_cast<double>(atom_index)});
        if (atom_index != 0) {
            bonds.emplace_back(atom_index - 1, atom_index, core::BondOrder::SINGLE);
        }
    }
    atoms.emplace_back(1, isolated_formal_charge, "isolated");
    positions.push_back(core::Position{.x = 30.0});

    return core::Molecule{std::move(atoms), std::move(bonds),
                          std::vector{core::Conformer{std::move(positions)}},
                          "extended-components"};
}

auto make_two_diatomic_components(const double separation) -> core::Molecule {
    return core::Molecule{
        std::vector{core::Atom{8, 0, "O1"}, core::Atom{1, 0, "H1"}, core::Atom{8, 0, "O2"},
                    core::Atom{1, 0, "H2"}},
        std::vector{core::Bond{0, 1, core::BondOrder::SINGLE},
                    core::Bond{2, 3, core::BondOrder::SINGLE}},
        std::vector{core::Conformer{{core::Position{.x = 0.0}, core::Position{.x = 1.0},
                                     core::Position{.x = separation},
                                     core::Position{.x = separation + 1.2}}}},
        "two-diatomic-components"};
}

[[nodiscard]] auto calculate_reduced(const core::Molecule& molecule,
                                     const std::string_view method_id,
                                     std::vector<parameters::ParameterSet> parameter_sets,
                                     const calculation::ExecutionSelectionKind mode,
                                     const double radius = 8.0) -> charges::AtomicCharges {
    const auto result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{molecule}},
        .parameter_sets = std::move(parameter_sets),
        .method_id = std::string{method_id},
        .execution_selection = calculation::ExecutionSelection{mode, radius}});
    REQUIRE(result.calculated());
    REQUIRE(result.charges.has_value());
    return result.charges->assignment(0).charges;
}

[[nodiscard]] auto calculate_reduced(const core::Molecule& molecule,
                                     const std::string_view method_id,
                                     const parameters::ParameterSet& parameter_set,
                                     const calculation::ExecutionSelectionKind mode,
                                     const double radius = 8.0) -> charges::AtomicCharges {
    return calculate_reduced(molecule, method_id, std::vector{parameter_set}, mode, radius);
}

[[nodiscard]] auto calculate_full(const core::Molecule& molecule, const std::string_view method_id,
                                  std::vector<parameters::ParameterSet> parameter_sets)
    -> charges::AtomicCharges {
    const auto result = calculate_application(
        calculation::AssessmentRequest{.molecules = core::MoleculeCollection{std::vector{molecule}},
                                       .parameter_sets = std::move(parameter_sets),
                                       .method_id = std::string{method_id}});
    REQUIRE(result.calculated());
    REQUIRE(result.charges.has_value());
    return result.charges->assignment(0).charges;
}

[[nodiscard]] auto calculate_full(const core::Molecule& molecule, const std::string_view method_id,
                                  const parameters::ParameterSet& parameter_set)
    -> charges::AtomicCharges {
    return calculate_full(molecule, method_id, std::vector{parameter_set});
}

struct ErrorMetrics {
    double rmsd;
    double mae;
    double maxabs;
};

[[nodiscard]] auto error_metrics(const charges::AtomicCharges& actual,
                                 const charges::AtomicCharges& reference) -> ErrorMetrics {
    REQUIRE(actual.size() == reference.size());
    auto squared_error = 0.0;
    auto absolute_error = 0.0;
    auto maxabs = 0.0;
    for (std::size_t atom_index = 0; atom_index < actual.size(); ++atom_index) {
        const auto difference = actual[atom_index] - reference[atom_index];
        const auto absolute_difference = std::abs(difference);
        squared_error += difference * difference;
        absolute_error += absolute_difference;
        maxabs = std::max(maxabs, absolute_difference);
    }
    const auto atom_count = static_cast<double>(actual.size());
    return {.rmsd = std::sqrt(squared_error / atom_count),
            .mae = absolute_error / atom_count,
            .maxabs = maxabs};
}

auto check_component_total(const charges::AtomicCharges& atomic_charges, const std::size_t begin,
                           const std::size_t end, const double expected) -> void {
    auto total = 0.0;
    for (auto atom_index = begin; atom_index < end; ++atom_index) {
        total += atomic_charges[atom_index];
    }
    CHECK(std::abs(total - expected) < 1.0e-10);
}

auto assert_reduced_matches_full(
    const std::string_view method_id,
    const std::vector<parameters::ParameterSet>& parameter_sets = {},
    core::Molecule molecule = chargefw::test::make_two_conformer_water()) -> void {
    const auto molecules = core::MoleculeCollection{std::vector{molecule, std::move(molecule)}};
    const auto full =
        calculate_application(calculation::AssessmentRequest{.molecules = molecules,
                                                             .parameter_sets = parameter_sets,
                                                             .method_id = std::string{method_id}});
    REQUIRE(full.calculated());

    for (const auto selection_kind : {calculation::ExecutionSelectionKind::cutoff,
                                      calculation::ExecutionSelectionKind::cover}) {
        const auto reduced = calculate_application(calculation::AssessmentRequest{
            .molecules = molecules,
            .parameter_sets = parameter_sets,
            .method_id = std::string{method_id},
            .execution_selection = calculation::ExecutionSelection{selection_kind, 8.0}});

        REQUIRE(reduced.calculated());
        REQUIRE(reduced.effective.has_value());
        CHECK(reduced.effective->execution_policy.mode() ==
              (selection_kind == calculation::ExecutionSelectionKind::cutoff
                   ? calculation::ExecutionMode::cutoff
                   : calculation::ExecutionMode::cover));
        REQUIRE(full.charges->size() == reduced.charges->size());

        for (std::size_t assignment_index = 0; assignment_index < full.charges->size();
             ++assignment_index) {
            const auto& full_charges = full.charges->assignment(assignment_index).charges;
            const auto& reduced_charges = reduced.charges->assignment(assignment_index).charges;
            REQUIRE(full_charges.size() == reduced_charges.size());
            for (std::size_t atom_index = 0; atom_index < full_charges.size(); ++atom_index) {
                CHECK(std::abs(full_charges[atom_index] - reduced_charges[atom_index]) < 1.0e-10);
            }
        }
    }
}

} // namespace

TEST_CASE("reduced execution validates inputs and mode selection",
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

    const auto no_conformer_collection = core::MoleculeCollection{std::vector{
        charged_molecule, core::Molecule{std::vector{core::Atom{1}}, {}, {}, "no-conformer"}}};
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
        try {
            calculate_without_conformer();
            CHECK(false);
        } catch (const std::invalid_argument& error) {
            CHECK(std::string_view{error.what()} ==
                  "selected method 'zero-fragment' requires coordinates, but molecule 2 has no "
                  "conformers");
        }
    }

    const auto cutoff = calculation::calculate(
        {.molecules = prepared,
         .selected = selected,
         .execution_policy = calculation::ExecutionPolicy{calculation::ExecutionMode::cutoff, 8.0},
         .max_threads = 2});
    CHECK(cutoff.charges.assignment(0).charges[0] == 0.5);
    CHECK(cutoff.charges.assignment(0).charges[1] == 0.5);

    const auto cover = calculation::calculate(
        {.molecules = prepared,
         .selected = selected,
         .execution_policy = calculation::ExecutionPolicy{calculation::ExecutionMode::cover, 8.0}});
    CHECK(cover.charges.assignment(0).charges[0] == 0.5);
    CHECK(cover.charges.assignment(0).charges[1] == 0.5);

    assert_reduced_matches_full("eem", {make_eem_parameters()});

    const auto automatic_cutoff = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {make_eem_parameters()},
        .method_id = "eem",
        .resource_policy = {.cutoff_atom_threshold = 2}});
    REQUIRE(automatic_cutoff.calculated());
    REQUIRE(automatic_cutoff.effective.has_value());
    CHECK(automatic_cutoff.effective->execution_policy.mode() ==
          calculation::ExecutionMode::cutoff);
    CHECK(automatic_cutoff.effective->execution_policy.radius() ==
          std::optional<double>{calculation::default_automatic_reduced_radius});

    const auto overridden_automatic_cutoff = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {make_eem_parameters()},
        .method_id = "eem",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::automatic, 8.0},
        .resource_policy = {.cutoff_atom_threshold = 2}});
    REQUIRE(overridden_automatic_cutoff.calculated());
    REQUIRE(overridden_automatic_cutoff.effective.has_value());
    CHECK(overridden_automatic_cutoff.effective->execution_policy.mode() ==
          calculation::ExecutionMode::cutoff);
    CHECK(overridden_automatic_cutoff.effective->execution_policy.radius() ==
          std::optional<double>{8.0});

    const auto automatic_cover = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {make_eem_parameters()},
        .method_id = "eem",
        .resource_policy = {.cutoff_atom_threshold = 2, .cover_atom_threshold = 2}});
    REQUIRE(automatic_cover.calculated());
    REQUIRE(automatic_cover.effective.has_value());
    CHECK(automatic_cover.effective->execution_policy.mode() == calculation::ExecutionMode::cover);

    const auto explicit_cutoff_above_cover_threshold =
        calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {make_eem_parameters()},
            .method_id = "eem",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff,
                                                calculation::minimum_reduced_radius},
            .resource_policy = {.cutoff_atom_threshold = 2, .cover_atom_threshold = 2}});
    REQUIRE(explicit_cutoff_above_cover_threshold.calculated());
    REQUIRE(explicit_cutoff_above_cover_threshold.effective.has_value());
    CHECK(explicit_cutoff_above_cover_threshold.effective->execution_policy.mode() ==
          calculation::ExecutionMode::cutoff);
    CHECK(explicit_cutoff_above_cover_threshold.effective->execution_issues.size() == 1);

    assert_reduced_matches_full("qeq", {make_qeq_parameters()});
    assert_reduced_matches_full("sfkeem", {make_sfkeem_parameters()});
    assert_reduced_matches_full("eqeq");
    assert_reduced_matches_full("eqeqc", {make_eqeqc_parameters()});
    assert_reduced_matches_full("abeem", {make_abeem_parameters()});
    assert_reduced_matches_full("sqe", {make_sqe_parameters("sqe", false)});
    assert_reduced_matches_full("sqeq0", {make_sqe_parameters("sqeq0", false)});
    assert_reduced_matches_full("sqeqp", {make_sqe_parameters("sqeqp", true)});
    assert_reduced_matches_full("sqe", {make_sqe_parameters("sqe", false, true)});
    assert_reduced_matches_full("sqeq0", {make_sqe_parameters("sqeq0", false, true)});
    assert_reduced_matches_full("sqeqp", {make_sqe_parameters("sqeqp", true, true)});
    assert_reduced_matches_full("sqeq0", {make_sqe_parameters("sqeq0", false)},
                                make_charged_water());
    assert_reduced_matches_full("sqeqp", {make_sqe_parameters("sqeqp", true)},
                                make_charged_water());
}

TEST_CASE("SQE reduced execution preserves original component charge budgets",
          "[calculation][reduced-execution][sqe]") {
    constexpr auto chain_atom_count = std::size_t{12};
    for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                            calculation::ExecutionSelectionKind::cover}) {
        CAPTURE(mode);

        const auto neutral = make_extended_components(0);
        const auto sqe = calculate_reduced(neutral, "sqe", make_sqe_parameters("sqe", false), mode);
        check_component_total(sqe, 0, chain_atom_count, 0.0);
        check_component_total(sqe, chain_atom_count, chain_atom_count + 1, 0.0);

        const auto charged = make_extended_components(1);
        const auto sqeq0 =
            calculate_reduced(charged, "sqeq0", make_sqe_parameters("sqeq0", false), mode);
        check_component_total(sqeq0, 0, chain_atom_count, 0.0);
        check_component_total(sqeq0, chain_atom_count, chain_atom_count + 1, 1.0);

        const auto sqeqp =
            calculate_reduced(charged, "sqeqp", make_sqe_parameters("sqeqp", true), mode);
        constexpr auto atom_count = chain_atom_count + 1;
        constexpr auto raw_total = -1.25;
        constexpr auto offset = (1.0 - raw_total) / static_cast<double>(atom_count);
        constexpr auto chain_reference_total = -1.5 + chain_atom_count * offset;
        constexpr auto isolated_reference_total = 0.25 + offset;
        check_component_total(sqeqp, 0, chain_atom_count, chain_reference_total);
        check_component_total(sqeqp, chain_atom_count, atom_count, isolated_reference_total);
    }
}

TEST_CASE("reduced approximation remains bounded across a truncated radius sweep",
          "[calculation][reduced-execution][accuracy]") {
    struct AccuracyCase {
        std::string_view method_id;
        core::Molecule molecule;
        std::vector<parameters::ParameterSet> parameter_sets;
        double max_rmsd;
        double max_mae;
        double max_maxabs;
    };
    const auto charged = make_extended_components(1);
    const auto cases =
        std::vector<AccuracyCase>{{.method_id = "abeem",
                                   .molecule = charged,
                                   .parameter_sets = {make_abeem_parameters()},
                                   .max_rmsd = 0.01,
                                   .max_mae = 0.01,
                                   .max_maxabs = 0.02},
                                  {.method_id = "eem",
                                   .molecule = charged,
                                   .parameter_sets = {make_eem_parameters()},
                                   .max_rmsd = 0.03,
                                   .max_mae = 0.02,
                                   .max_maxabs = 0.08},
                                  {.method_id = "eqeq",
                                   .molecule = charged,
                                   .parameter_sets = {},
                                   .max_rmsd = 0.05,
                                   .max_mae = 0.03,
                                   .max_maxabs = 0.15},
                                  {.method_id = "eqeqc",
                                   .molecule = charged,
                                   .parameter_sets = {make_eqeqc_parameters()},
                                   .max_rmsd = 0.05,
                                   .max_mae = 0.03,
                                   .max_maxabs = 0.15},
                                  {.method_id = "qeq",
                                   .molecule = charged,
                                   .parameter_sets = {make_qeq_parameters()},
                                   .max_rmsd = 0.1,
                                   .max_mae = 0.05,
                                   .max_maxabs = 0.3},
                                  {.method_id = "sfkeem",
                                   .molecule = charged,
                                   .parameter_sets = {make_sfkeem_parameters()},
                                   .max_rmsd = 0.04,
                                   .max_mae = 0.02,
                                   .max_maxabs = 0.1},
                                  {.method_id = "sqe",
                                   .molecule = make_extended_components(0),
                                   .parameter_sets = {make_sqe_parameters("sqe", false)},
                                   .max_rmsd = 0.01,
                                   .max_mae = 0.01,
                                   .max_maxabs = 0.02},
                                  {.method_id = "sqeq0",
                                   .molecule = charged,
                                   .parameter_sets = {make_sqe_parameters("sqeq0", false)},
                                   .max_rmsd = 0.01,
                                   .max_mae = 0.01,
                                   .max_maxabs = 0.02},
                                  {.method_id = "sqeqp",
                                   .molecule = charged,
                                   .parameter_sets = {make_sqe_parameters("sqeqp", true)},
                                   .max_rmsd = 0.01,
                                   .max_mae = 0.01,
                                   .max_maxabs = 0.02}};

    auto aggregate_squared_error = 0.0;
    auto aggregate_absolute_error = 0.0;
    auto aggregate_maxabs = 0.0;
    auto aggregate_atom_count = std::size_t{0};

    for (const auto& test_case : cases) {
        const features::PreparedMolecule prepared{test_case.molecule};
        const features::ConformerFeatures geometry{test_case.molecule, 0};
        const features::SpatialFragmentBuilder builder{prepared, geometry};
        CHECK(builder.build(0, 8.0).molecule().atom_count() < test_case.molecule.atom_count());

        const auto full =
            calculate_full(test_case.molecule, test_case.method_id, test_case.parameter_sets);
        for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                                calculation::ExecutionSelectionKind::cover}) {
            for (const auto radius : {8.0, 12.0, 16.0}) {
                const auto reduced = calculate_reduced(test_case.molecule, test_case.method_id,
                                                       test_case.parameter_sets, mode, radius);
                const auto metrics = error_metrics(reduced, full);
                aggregate_squared_error +=
                    metrics.rmsd * metrics.rmsd * static_cast<double>(reduced.size());
                aggregate_absolute_error += metrics.mae * static_cast<double>(reduced.size());
                aggregate_maxabs = std::max(aggregate_maxabs, metrics.maxabs);
                aggregate_atom_count += reduced.size();
                CAPTURE(test_case.method_id, mode, radius, metrics.rmsd, metrics.mae,
                        metrics.maxabs);
                CHECK(std::isfinite(metrics.rmsd));
                CHECK(std::isfinite(metrics.mae));
                CHECK(std::isfinite(metrics.maxabs));
                CHECK(metrics.rmsd < test_case.max_rmsd);
                CHECK(metrics.mae < test_case.max_mae);
                CHECK(metrics.maxabs < test_case.max_maxabs);
            }
        }
    }

    const auto aggregate_count = static_cast<double>(aggregate_atom_count);
    const auto aggregate_rmsd = std::sqrt(aggregate_squared_error / aggregate_count);
    const auto aggregate_mae = aggregate_absolute_error / aggregate_count;
    CAPTURE(aggregate_rmsd, aggregate_mae, aggregate_maxabs, aggregate_atom_count);
    CHECK(aggregate_rmsd < 0.05);
    CHECK(aggregate_mae < 0.03);
    CHECK(aggregate_maxabs < 0.3);
}

TEST_CASE("SQE component totals coexist with intermolecular polarization",
          "[calculation][reduced-execution][sqe]") {
    const auto parameters = make_sqe_parameters("sqeq0", false);
    const auto near = make_two_diatomic_components(4.0);
    const auto far = make_two_diatomic_components(30.0);

    const auto full_near = calculate_full(near, "sqeq0", parameters);
    const auto full_far = calculate_full(far, "sqeq0", parameters);
    check_component_total(full_near, 0, 2, 0.0);
    check_component_total(full_near, 2, 4, 0.0);
    CHECK(std::abs(full_near[0] - full_far[0]) > 1.0e-6);

    for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                            calculation::ExecutionSelectionKind::cover}) {
        const auto reduced_near = calculate_reduced(near, "sqeq0", parameters, mode);
        const auto reduced_far = calculate_reduced(far, "sqeq0", parameters, mode);
        CAPTURE(mode);
        check_component_total(reduced_near, 0, 2, 0.0);
        check_component_total(reduced_near, 2, 4, 0.0);
        CHECK(std::abs(reduced_near[0] - reduced_far[0]) > 1.0e-6);
    }
}

TEST_CASE("reduced solver failures retain method and target context",
          "[calculation][reduced-execution]") {
    for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                            calculation::ExecutionSelectionKind::cover}) {
        const auto result = calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {make_invalid_qeq_parameters()},
            .method_id = "qeq",
            .parameter_set_id = "invalid-qeq",
            .execution_selection =
                calculation::ExecutionSelection{mode, calculation::minimum_reduced_radius}});

        CHECK(result.status == calculation::ExecutionStatus::numerical_failure);
        CHECK_FALSE(result.calculated());
        REQUIRE(result.failure_message.has_value());
        const auto message = std::string_view{*result.failure_message};
        CHECK(message.contains("method 'qeq', molecule 1 ('water'), conformer 1 failed:"));
        CHECK(message.contains(mode == calculation::ExecutionSelectionKind::cutoff
                                   ? "cutoff fragment around source atom 1 failed:"
                                   : "cover fragment around source atom 1 failed:"));
        CHECK_FALSE(message.contains("center atom"));
        CHECK_FALSE(message.contains("pivot atom"));
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
                 calculation::ExecutionPolicy{mode, calculation::minimum_reduced_radius},
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
