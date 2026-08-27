#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <algorithm>
#include <chargefw/calculation/calculation.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <concepts>
#include <limits>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

template <typename T>
concept HasPublicParameterSets = requires(T& value) { value.parameter_sets; };

template <typename T>
concept HasPublicSelectedCandidate = requires(T& value) { value.selected; };

[[nodiscard]] auto calculate_application(const calculation::AssessmentRequest& request)
    -> calculation::ExecutionResult {
    const auto max_threads = request.resource_policy.max_threads;
    auto assessment = calculation::assess(request);
    return calculation::calculate(std::move(assessment), max_threads);
}

auto make_parameter_set(std::string id, std::string method_id, const std::uint16_t priority)
    -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{.id = std::move(id),
                                                   .method_id = std::move(method_id),
                                                   .name = "Test parameters",
                                                   .priority = priority},
        {},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, chargefw::parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, chargefw::parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "value", .value = 1.0}}}}}};
}

auto make_singular_eem_parameter_set() -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{
            .id = "singular-eem", .method_id = "eem", .name = "Singular EEM parameters"},
        chargefw::parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, chargefw::parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, chargefw::parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 1.0}}}}}};
}

auto make_singular_eem_molecule() -> core::Molecule {
    return core::Molecule{
        std::vector{core::Atom{1}, core::Atom{8}},
        {},
        {core::Conformer{{core::Position{0.0, 0.0, 0.0}, core::Position{1.0, 0.0, 0.0}}}},
        "singular-eem"};
}

auto make_duplicate_parameter_request(const std::string_view first_method_id,
                                      const std::string_view second_method_id,
                                      const bool explicit_selection)
    -> calculation::AssessmentRequest {
    auto request = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {make_parameter_set("duplicate", std::string{first_method_id}, 0),
                           make_parameter_set("duplicate", std::string{second_method_id}, 0)}};
    if (explicit_selection) {
        request.method_id = "qeq";
        request.parameter_set_id = "duplicate";
    }
    return request;
}

auto make_double_bonded_carbons() -> core::Molecule {
    return core::Molecule{std::vector{core::Atom{6}, core::Atom{6}},
                          std::vector{core::Bond{0, 1, core::BondOrder::DOUBLE}},
                          {},
                          "double-bonded-carbons"};
}

auto make_permissive_parameter_set() -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{.id = "permissive-parameters",
                                                   .method_id = "permissive-parameterized",
                                                   .name = "Permissive parameters"},
        {},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  6, chargefw::parameters::AtomParameterClassificationKind::HIGHEST_BOND_ORDER,
                  "1"),
              .parameters = {{.name = "value", .value = 3.0}}}}}};
}

auto make_permissive_peoe_parameter_set() -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{.id = "permissive-peoe-parameters",
                                                   .method_id = "peoe",
                                                   .name = "Permissive PEOE parameters"},
        chargefw::parameters::CommonParameters{{{.name = "dampH", .value = 1.0}}},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  6, chargefw::parameters::AtomParameterClassificationKind::HIGHEST_BOND_ORDER,
                  "1"),
              .parameters = {{.name = "A", .value = 1.0},
                             {.name = "B", .value = 1.0},
                             {.name = "C", .value = 1.0}}}}}};
}

} // namespace

static_assert(!HasPublicParameterSets<calculation::AssessmentResult>);
static_assert(!HasPublicSelectedCandidate<calculation::AssessmentResult>);
static_assert(std::is_move_constructible_v<calculation::AssessmentResult>);
static_assert(!std::is_move_assignable_v<calculation::AssessmentResult>);

TEST_CASE("calculation facade reports singular solver failures with target context",
          "[calculation][calculation]") {
    const auto result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_singular_eem_molecule()}},
        .parameter_sets = {make_singular_eem_parameter_set()},
        .method_id = "eem",
        .parameter_set_id = "singular-eem"});

    CHECK(result.status == calculation::ExecutionStatus::numerical_failure);
    CHECK(!result.charges.has_value());
    REQUIRE(result.effective.has_value());
    CHECK(result.effective->method_id == "eem");
    CHECK(result.effective->parameter_set_id == std::optional<std::string>{"singular-eem"});
    REQUIRE(result.failure_message.has_value());
    CHECK(result.failure_message->find("method 'eem', molecule 1 ('singular-eem'), conformer 1") !=
          std::string::npos);
}

TEST_CASE("assessment preserves owned selection state and validates method options",
          "[calculation][calculation]") {

    // The returned report owns its IDs and does not retain parameter-set pointers from the consumed
    // assessment. Moving the assessment before execution must preserve the indexed selected plan.
    const auto owned_parameterized_result = []() -> calculation::ExecutionResult {
        auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
            .parameter_sets = {make_permissive_peoe_parameter_set()},
            .method_id = "peoe",
            .parameter_set_id = "permissive-peoe-parameters",
            .classification_options = {.permissive_types = true}});
        auto relocated_assessment = std::move(assessment);
        return calculation::calculate(std::move(relocated_assessment));
    }();
    REQUIRE(owned_parameterized_result.calculated());
    REQUIRE(owned_parameterized_result.applicability.applicable.size() == 1);
    CHECK(owned_parameterized_result.applicability.applicable[0].method_id == "peoe");
    CHECK(owned_parameterized_result.applicability.applicable[0].parameter_set_id ==
          std::optional<std::string>{"permissive-peoe-parameters"});
    CHECK(owned_parameterized_result.applicability.selected_candidate_index ==
          std::optional<std::size_t>{0});

    // Rvalue assessment transfers the request's heavy inputs while preserving its lightweight
    // selection data until assessment has completed.
    auto consuming_request = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water(),
                                                          chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "formal",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}};
    auto consuming_assessment = calculation::assess(std::move(consuming_request));
    REQUIRE(consuming_assessment.executable());
    CHECK(consuming_assessment.applicability().selected_candidate_index ==
          std::optional<std::size_t>{0});
    const auto consuming_result = calculation::calculate(std::move(consuming_assessment));
    REQUIRE(consuming_result.calculated());
    REQUIRE(consuming_result.charges->size() == 2);
    CHECK(consuming_result.charges->assignment(0).target.molecule_index == 0);
    CHECK(consuming_result.charges->assignment(1).target.molecule_index == 1);

    auto peoe_options = methods::MethodOptions{};
    peoe_options.set("iters", 1);
    const auto configured_peoe_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
        .parameter_sets = {make_permissive_peoe_parameter_set()},
        .method_id = "peoe",
        .parameter_set_id = "permissive-peoe-parameters",
        .method_options = {{"peoe", peoe_options}},
        .classification_options = {.permissive_types = true}});
    REQUIRE(configured_peoe_result.calculated());
    REQUIRE(configured_peoe_result.effective.has_value());
    CHECK(configured_peoe_result.effective->method_options.get<int>("iters") == 1);

    auto invalid_peoe_options = methods::MethodOptions{};
    invalid_peoe_options.set("iters", std::string{"one"});
    const auto calculate_with_invalid_peoe_options = [&] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
            .parameter_sets = {make_permissive_peoe_parameter_set()},
            .method_id = "peoe",
            .parameter_set_id = "permissive-peoe-parameters",
            .method_options = {{"peoe", invalid_peoe_options}},
            .classification_options = {.permissive_types = true}}));
    };
    CHECK_THROWS_AS(calculate_with_invalid_peoe_options(), std::invalid_argument);

    const auto rejected_explicit_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "smpqeq"});
    CHECK(!rejected_explicit_assessment.executable());
    REQUIRE(rejected_explicit_assessment.applicability().rejected.size() == 1);
    CHECK(rejected_explicit_assessment.applicability().rejected[0].method_id == "smpqeq");
    CHECK(!rejected_explicit_assessment.applicability().rejected[0].issues.empty());

    const auto rejected_parameter_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
        .parameter_sets = {make_permissive_peoe_parameter_set()},
        .method_id = "peoe",
        .parameter_set_id = "permissive-peoe-parameters"});
    CHECK(!rejected_parameter_assessment.executable());
    REQUIRE(rejected_parameter_assessment.applicability().rejected.size() == 1);
    CHECK(rejected_parameter_assessment.applicability().rejected[0].method_id == "peoe");
    CHECK(rejected_parameter_assessment.applicability().rejected[0].parameter_set_id ==
          std::optional<std::string>{"permissive-peoe-parameters"});

    auto automatic_rejected_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}}});
    const auto automatic_rejected_result =
        calculation::calculate(std::move(automatic_rejected_assessment));
    const auto rejected_smpqeq =
        std::ranges::find_if(automatic_rejected_result.applicability.rejected,
                             [](const calculation::RejectedCandidateReport& candidate) {
                                 return candidate.method_id == "smpqeq";
                             });
    REQUIRE(rejected_smpqeq != automatic_rejected_result.applicability.rejected.end());
    CHECK(!rejected_smpqeq->issues.empty());
}

TEST_CASE("calculation facade applies execution policy and rejects invalid plans",
          "[calculation][calculation]") {

    const auto application_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "formal",
        .parameter_set_id = std::nullopt,
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});

    REQUIRE(application_result.charges.has_value());
    CHECK(application_result.charges->method_id() == std::string_view{"formal"});
    CHECK(application_result.charges->size() == 1);
    REQUIRE(application_result.effective.has_value());
    CHECK(application_result.effective->method_id == "formal");
    CHECK_FALSE(application_result.effective->parameter_set_id.has_value());
    CHECK(application_result.effective->execution_policy.mode() ==
          calculation::ExecutionMode::full);
    CHECK(application_result.effective->execution_issues.empty());

    auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "formal",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});
    REQUIRE(assessment.executable());
    REQUIRE(assessment.applicability().applicable.size() == 1);
    REQUIRE(assessment.execution_policy().has_value());
    CHECK(assessment.applicability().selected_candidate_index == std::optional<std::size_t>{0});
    CHECK(assessment.applicability().applicable[0].method_id == "formal");
    CHECK(assessment.execution_policy()->mode() == calculation::ExecutionMode::full);

    const auto assessed_result = calculation::calculate(std::move(assessment), 1);
    REQUIRE(assessed_result.calculated());
    CHECK(assessed_result.charges->method_id() == std::string_view{"formal"});
    CHECK(assessed_result.metrics.applicability_seconds >= 0.0);

    const auto automatic_fallback_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .resource_policy = {.cutoff_atom_threshold = 2}});
    REQUIRE(automatic_fallback_result.calculated());
    CHECK(automatic_fallback_result.charges->method_id() == std::string_view{"eqeq"});
    REQUIRE(automatic_fallback_result.effective.has_value());
    CHECK(automatic_fallback_result.effective->method_id == "eqeq");
    CHECK(automatic_fallback_result.effective->execution_policy.mode() ==
          calculation::ExecutionMode::cutoff);
    CHECK(automatic_fallback_result.effective->execution_policy.radius() ==
          std::optional<double>{calculation::default_automatic_reduced_radius});
    CHECK(automatic_fallback_result.effective->execution_issues.empty());

    const auto automatic_mgc_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "mgc",
        .resource_policy = {.cutoff_atom_threshold = 2}});
    CHECK(automatic_mgc_result.status == calculation::ExecutionStatus::no_executable_plan);
    CHECK_FALSE(automatic_mgc_result.calculated());
    CHECK_FALSE(automatic_mgc_result.effective.has_value());

    const auto explicit_full_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "mgc",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
        .resource_policy = {.cutoff_atom_threshold = 2}});
    REQUIRE(explicit_full_result.calculated());
    CHECK(explicit_full_result.charges->method_id() == std::string_view{"mgc"});
    REQUIRE(explicit_full_result.effective.has_value());
    CHECK(explicit_full_result.effective->execution_policy.mode() ==
          calculation::ExecutionMode::full);
    REQUIRE(explicit_full_result.effective->execution_issues.size() == 1);
    CHECK(explicit_full_result.effective->execution_issues[0].kind ==
          methods::ExecutionIssueKind::resource_threshold_exceeded);

    const auto unlimited_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "mgc",
        .resource_policy = {.cutoff_atom_threshold = std::nullopt,
                            .cover_atom_threshold = std::nullopt}});
    REQUIRE(unlimited_result.calculated());
    CHECK(unlimited_result.charges->method_id() == std::string_view{"mgc"});
    REQUIRE(unlimited_result.effective.has_value());
    CHECK(unlimited_result.effective->execution_issues.empty());

    const auto unsupported_cutoff_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "mgc",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff, 8.0}});
    CHECK(unsupported_cutoff_result.status == calculation::ExecutionStatus::no_executable_plan);
    CHECK_FALSE(unsupported_cutoff_result.calculated());

    const auto calculate_missing_method = [] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "missing",
            .parameter_set_id = std::nullopt}));
    };
    CHECK_THROWS_AS(calculate_missing_method(), std::invalid_argument);

    const auto calculate_missing_parameter_set = [] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .parameter_set_id = "missing"}));
    };
    CHECK_THROWS_AS(calculate_missing_parameter_set(), std::invalid_argument);

    const auto assess_missing_cutoff_threshold = [] -> void {
        static_cast<void>(calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .resource_policy = {.cutoff_atom_threshold = std::nullopt,
                                .cover_atom_threshold = 10}}));
    };
    CHECK_THROWS_AS(assess_missing_cutoff_threshold(), std::invalid_argument);
    const auto assess_inconsistent_thresholds = [] -> void {
        static_cast<void>(calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .resource_policy = {.cutoff_atom_threshold = 20, .cover_atom_threshold = 10}}));
    };
    CHECK_THROWS_AS(assess_inconsistent_thresholds(), std::invalid_argument);
}

TEST_CASE("calculation preserves empty-input cardinality and assessment ownership",
          "[calculation][calculation]") {

    // Empty collections and empty targets retain their distinct cardinalities: no collection
    // entries produce no assignments, while an empty molecule still produces one source assignment.
    const auto empty_collection_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector<core::Molecule>{}},
        .method_id = "formal"});
    REQUIRE(empty_collection_result.calculated());
    CHECK(empty_collection_result.charges->empty());

    const auto empty_target_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{core::Molecule{{}, {}, {}, "empty"}}},
        .method_id = "formal"});
    REQUIRE(empty_target_result.calculated());
    REQUIRE(empty_target_result.charges->size() == 1);
    CHECK(empty_target_result.charges->assignment(0).charges.empty());

    // AssessmentRequest rejects duplicate parameter-set IDs before filtering or applicability. This
    // remains true when the duplicates target the same or different methods and for both overloads.
    for (const auto explicit_selection : {false, true}) {
        for (const auto [first_method_id, second_method_id] :
             {std::pair{"qeq", "qeq"}, std::pair{"qeq", "eem"}}) {
            const auto lvalue_request = make_duplicate_parameter_request(
                first_method_id, second_method_id, explicit_selection);
            const auto assess_duplicate_lvalue_request = [&lvalue_request] -> void {
                static_cast<void>(calculation::assess(lvalue_request));
            };
            CHECK_THROWS_AS(assess_duplicate_lvalue_request(), std::invalid_argument);

            auto rvalue_request = make_duplicate_parameter_request(
                first_method_id, second_method_id, explicit_selection);
            const auto assess_duplicate_rvalue_request = [&rvalue_request] -> void {
                static_cast<void>(calculation::assess(std::move(rvalue_request)));
            };
            CHECK_THROWS_AS(assess_duplicate_rvalue_request(), std::invalid_argument);
        }
    }

    // Both assessment overloads retain the same selection state. The lvalue path copies its inputs,
    // so execution remains valid after the original request has been destroyed.
    auto lvalue_assessment = []() -> calculation::AssessmentResult {
        auto request = calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .method_id = "formal",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}};
        return calculation::assess(request);
    }();
    auto rvalue_request = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "formal",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}};
    auto rvalue_assessment = calculation::assess(std::move(rvalue_request));
    REQUIRE(lvalue_assessment.execution_policy().has_value());
    REQUIRE(rvalue_assessment.execution_policy().has_value());
    CHECK(lvalue_assessment.applicability().selected_candidate_index ==
          rvalue_assessment.applicability().selected_candidate_index);
    CHECK(lvalue_assessment.execution_policy()->mode() ==
          rvalue_assessment.execution_policy()->mode());
    const auto lvalue_result = calculation::calculate(std::move(lvalue_assessment));
    const auto rvalue_result = calculation::calculate(std::move(rvalue_assessment));
    REQUIRE(lvalue_result.calculated());
    REQUIRE(rvalue_result.calculated());
    REQUIRE(lvalue_result.charges->size() == 1);
    REQUIRE(rvalue_result.charges->size() == 1);
    CHECK(std::ranges::equal(lvalue_result.charges->assignment(0).charges.values(),
                             rvalue_result.charges->assignment(0).charges.values()));
}

TEST_CASE("parallel calculation materializes assignments in source order",
          "[calculation][calculation]") {

    // Parallel execution may emit progress out of order, but materialized assignments always retain
    // the source molecule/conformer order.
    const auto parallel_ordered_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{
            chargefw::test::make_two_conformer_water(), chargefw::test::make_water()}},
        .method_id = "eqeq",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
        .resource_policy = {.max_threads = 2}});
    REQUIRE(parallel_ordered_result.calculated());
    REQUIRE(parallel_ordered_result.charges->size() == 3);
    CHECK(parallel_ordered_result.charges->assignment(0).target.molecule_index == 0);
    CHECK(parallel_ordered_result.charges->assignment(0).target.conformer_index ==
          std::optional<std::size_t>{0});
    CHECK(parallel_ordered_result.charges->assignment(1).target.molecule_index == 0);
    CHECK(parallel_ordered_result.charges->assignment(1).target.conformer_index ==
          std::optional<std::size_t>{1});
    CHECK(parallel_ordered_result.charges->assignment(2).target.molecule_index == 1);
    CHECK(parallel_ordered_result.charges->assignment(2).target.conformer_index ==
          std::optional<std::size_t>{0});
}
