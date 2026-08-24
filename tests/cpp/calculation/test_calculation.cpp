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

class FixedChargeMethod : public methods::Method {
  public:
    FixedChargeMethod(std::string_view id, const int priority, const double value)
        : metadata_{.id = id,
                    .name = id,
                    .full_name = id,
                    .publication = std::nullopt,
                    .priority = priority},
          value_{value} {}

    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        return metadata_;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        return {};
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> charges::AtomicCharges override {
        return charges::AtomicCharges{std::vector<double>(input.molecule().atom_count(), value_)};
    }

  private:
    methods::MethodMetadata metadata_;
    double value_;
};

class ParameterizedFixedChargeMethod final : public FixedChargeMethod {
  public:
    ParameterizedFixedChargeMethod(std::string_view id, const int priority, const double value)
        : FixedChargeMethod{id, priority, value} {}

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        auto result = methods::MethodRequirements{};
        result.atom_parameters = {"value"};
        return result;
    }
};

class ParameterValueMethod final : public FixedChargeMethod {
  public:
    ParameterValueMethod() : FixedChargeMethod{"permissive-parameterized", 0, 0.0} {}

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        auto result = methods::MethodRequirements{};
        result.atom_parameters = {"value"};
        return result;
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> charges::AtomicCharges override {
        const auto values = input.parameters().atom("value");
        auto result = std::vector<double>{};
        result.reserve(input.molecule().atom_count());
        for (std::size_t atom_index = 0; atom_index < input.molecule().atom_count(); ++atom_index) {
            result.push_back(values[atom_index]);
        }
        return charges::AtomicCharges{std::move(result)};
    }
};

auto make_prepared_water() -> features::PreparedMoleculeCollection {
    static const core::MoleculeCollection collection{std::vector{chargefw::test::make_water()}};
    return features::PreparedMoleculeCollection{collection};
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

struct AutomaticCalculationResult {
    std::optional<charges::ChargeSet> charges;
    methods::ApplicabilityResult applicability;

    [[nodiscard]] auto calculated() const noexcept -> bool {
        return charges.has_value();
    }
};

auto calculate_automatically(
    const features::PreparedMoleculeCollection& molecules,
    std::span<const methods::Method* const> candidate_methods,
    std::span<const chargefw::parameters::ParameterSet> parameter_sets,
    const chargefw::parameters::ClassificationOptions classification_options = {})
    -> AutomaticCalculationResult {
    auto applicability =
        methods::find_applicable_methods({.molecules = molecules,
                                          .methods = candidate_methods,
                                          .parameter_sets = parameter_sets,
                                          .classification_options = classification_options});
    const auto* selected = calculation::select_applicable_method(applicability);

    if (selected == nullptr) {
        return {.charges = std::nullopt, .applicability = std::move(applicability)};
    }

    auto result = calculation::calculate({.molecules = molecules, .selected = *selected});
    return {.charges = std::move(result.charges), .applicability = std::move(applicability)};
}

} // namespace

static_assert(!HasPublicParameterSets<calculation::AssessmentResult>);
static_assert(!HasPublicSelectedCandidate<calculation::AssessmentResult>);
static_assert(std::is_move_constructible_v<calculation::AssessmentResult>);
static_assert(!std::is_move_assignable_v<calculation::AssessmentResult>);

TEST_CASE("calculation selects the highest-priority full execution plan",
          "[calculation][calculation]") {
    const auto prepared = make_prepared_water();
    const FixedChargeMethod higher_priority{"higher", 10, 10.0};
    const FixedChargeMethod lower_priority{"lower", 1, 1.0};
    const std::vector<const methods::Method*> methods{&lower_priority, &higher_priority};
    const std::vector<chargefw::parameters::ParameterSet> parameters;

    const auto applicability = methods::find_applicable_methods(
        {.molecules = prepared, .methods = methods, .parameter_sets = parameters});
    const auto* selected = calculation::select_applicable_method(applicability);
    CHECK(selected != nullptr);
    const auto plan =
        calculation::select_execution_plan(applicability, calculation::ExecutionSelection{});
    CHECK(plan.has_value());
    CHECK(plan->selected == selected);
    CHECK(plan->policy.mode() == calculation::ExecutionMode::full);
    CHECK(plan->issues.empty());
    const auto result = calculation::calculate({.molecules = prepared, .selected = *selected});
    CHECK(result.charges.method_id() == std::string_view{"higher"});
    CHECK(result.charges.size() == 1);
    CHECK(result.charges.assignment(0).charges[0] == 10.0);
    CHECK(applicability.applicable.size() == 2);

    const auto parallel_collection = core::MoleculeCollection{
        std::vector{chargefw::test::make_water(), chargefw::test::make_water()}};
    const features::PreparedMoleculeCollection parallel_prepared{parallel_collection};
    const auto parallel_result =
        calculation::calculate({.molecules = parallel_prepared,
                                .selected = *selected,
                                .execution_policy = calculation::ExecutionPolicy{},
                                .max_threads = 2});
    CHECK(parallel_result.charges.size() == 2);
    CHECK(parallel_result.charges.assignment(0).target.molecule_index == 0);
    CHECK(parallel_result.charges.assignment(1).target.molecule_index == 1);
    CHECK(parallel_result.charges.assignment(0).charges[0] == 10.0);
    CHECK(parallel_result.charges.assignment(1).charges[0] == 10.0);

    const auto throw_fn_0 = [&] -> void {
        static_cast<void>(
            calculation::calculate({.molecules = prepared,
                                    .selected = *selected,
                                    .execution_policy = calculation::ExecutionPolicy{},
                                    .max_threads = std::numeric_limits<std::size_t>::max()}));
    };
    CHECK_THROWS_AS(throw_fn_0(), std::invalid_argument);
}

TEST_CASE("automatic calculation ranks candidates and supports permissive classification",
          "[calculation][calculation]") {
    const auto prepared = make_prepared_water();
    const std::vector<chargefw::parameters::ParameterSet> parameters;

    const FixedChargeMethod alpha{"alpha", 1, 2.0};
    const FixedChargeMethod beta{"beta", 1, 3.0};
    const std::vector<const methods::Method*> tied_methods{&beta, &alpha};

    const auto tied_result = calculate_automatically(prepared, tied_methods, parameters);

    REQUIRE(tied_result.charges.has_value());
    CHECK(tied_result.charges->method_id() == std::string_view{"alpha"});
    CHECK(tied_result.charges->assignment(0).charges[0] == 2.0);

    const ParameterizedFixedChargeMethod parameterized{"parameterized", 0, 4.0};
    const std::vector<const methods::Method*> parameterized_methods{&parameterized};
    const std::vector parameter_sets{make_parameter_set("alpha", "parameterized", 1),
                                     make_parameter_set("zeta", "parameterized", 10)};

    const auto parameterized_result =
        calculate_automatically(prepared, parameterized_methods, parameter_sets);

    REQUIRE(parameterized_result.charges.has_value());
    CHECK(parameterized_result.charges->method_id() == std::string_view{"parameterized"});
    CHECK(parameterized_result.charges->parameter_set_id() == std::string_view{"zeta"});

    const std::vector<const methods::Method*> no_methods;
    const auto no_result = calculate_automatically(prepared, no_methods, parameters);

    CHECK(!no_result.calculated());
    CHECK(no_result.applicability.empty());

    const core::MoleculeCollection double_bonded_collection{
        std::vector{make_double_bonded_carbons()}};
    const features::PreparedMoleculeCollection double_bonded_prepared{double_bonded_collection};
    const ParameterValueMethod permissive_method;
    const std::vector<const methods::Method*> permissive_methods{&permissive_method};
    const std::vector permissive_parameter_sets{make_permissive_parameter_set()};

    const auto strict_permissive_result = calculate_automatically(
        double_bonded_prepared, permissive_methods, permissive_parameter_sets);
    CHECK(!strict_permissive_result.calculated());

    const auto permissive_calculation_result =
        calculate_automatically(double_bonded_prepared, permissive_methods,
                                permissive_parameter_sets, {.permissive_types = true});
    CHECK(permissive_calculation_result.calculated());
    CHECK(permissive_calculation_result.applicability.applicable.size() == 1);
    CHECK(permissive_calculation_result.charges->assignment(0).charges[0] == 3.0);
    CHECK(permissive_calculation_result.charges->assignment(0).charges[1] == 3.0);

    const auto throw_fn_1 = [] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
            .parameter_sets = {make_permissive_peoe_parameter_set()},
            .method_id = "peoe",
            .parameter_set_id = "permissive-peoe-parameters"}));
    };
    CHECK_THROWS_AS(throw_fn_1(), std::invalid_argument);

    const auto permissive_application_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
        .parameter_sets = {make_permissive_peoe_parameter_set()},
        .method_id = "peoe",
        .parameter_set_id = "permissive-peoe-parameters",
        .classification_options = {.permissive_types = true}});
    CHECK(permissive_application_result.calculated());
    CHECK(permissive_application_result.charges->method_id() == std::string_view{"peoe"});
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
    CHECK(owned_parameterized_result.calculated());
    CHECK(owned_parameterized_result.applicability.applicable.size() == 1);
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
    CHECK(consuming_assessment.executable());
    CHECK(consuming_assessment.applicability().selected_candidate_index ==
          std::optional<std::size_t>{0});
    const auto consuming_result = calculation::calculate(std::move(consuming_assessment));
    CHECK(consuming_result.calculated());
    CHECK(consuming_result.charges->size() == 2);
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
    CHECK(configured_peoe_result.calculated());
    CHECK(configured_peoe_result.effective_method_options.has_value());
    CHECK(configured_peoe_result.effective_method_options->get<int>("iters") == 1);

    auto invalid_peoe_options = methods::MethodOptions{};
    invalid_peoe_options.set("iters", std::string{"one"});
    const auto throw_fn_2 = [&] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
            .parameter_sets = {make_permissive_peoe_parameter_set()},
            .method_id = "peoe",
            .parameter_set_id = "permissive-peoe-parameters",
            .method_options = {{"peoe", invalid_peoe_options}},
            .classification_options = {.permissive_types = true}}));
    };
    CHECK_THROWS_AS(throw_fn_2(), std::invalid_argument);

    const auto rejected_explicit_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "smpqeq"});
    CHECK(!rejected_explicit_assessment.executable());
    CHECK(rejected_explicit_assessment.applicability().rejected.size() == 1);
    CHECK(rejected_explicit_assessment.applicability().rejected[0].method_id == "smpqeq");
    CHECK(!rejected_explicit_assessment.applicability().rejected[0].issues.empty());

    const auto rejected_parameter_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
        .parameter_sets = {make_permissive_peoe_parameter_set()},
        .method_id = "peoe",
        .parameter_set_id = "permissive-peoe-parameters"});
    CHECK(!rejected_parameter_assessment.executable());
    CHECK(rejected_parameter_assessment.applicability().rejected.size() == 1);
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
    CHECK(rejected_smpqeq != automatic_rejected_result.applicability.rejected.end());
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
    CHECK(application_result.execution_policy.has_value());
    CHECK(application_result.execution_policy->mode() == calculation::ExecutionMode::full);
    CHECK(application_result.execution_issues.empty());

    auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "formal",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});
    CHECK(assessment.executable());
    CHECK(assessment.applicability().applicable.size() == 1);
    CHECK(assessment.applicability().selected_candidate_index == std::optional<std::size_t>{0});
    CHECK(assessment.applicability().applicable[0].method_id == "formal");
    CHECK(assessment.execution_policy()->mode() == calculation::ExecutionMode::full);

    const auto assessed_result = calculation::calculate(std::move(assessment), 1);
    CHECK(assessed_result.calculated());
    CHECK(assessed_result.charges->method_id() == std::string_view{"formal"});
    CHECK(assessed_result.metrics.applicability_seconds >= 0.0);

    const auto automatic_fallback_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .resource_policy = {.cutoff_atom_threshold = 2}});
    CHECK(automatic_fallback_result.calculated());
    CHECK(automatic_fallback_result.charges->method_id() == std::string_view{"eqeq"});
    CHECK(automatic_fallback_result.execution_policy.has_value());
    CHECK(automatic_fallback_result.execution_policy->mode() == calculation::ExecutionMode::cutoff);
    CHECK(automatic_fallback_result.execution_policy->radius() ==
          std::optional<double>{calculation::default_automatic_reduced_radius});
    CHECK(automatic_fallback_result.execution_issues.empty());

    const auto throw_fn_3 = [] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "mgc",
            .resource_policy = {.cutoff_atom_threshold = 2}}));
    };
    CHECK_THROWS_AS(throw_fn_3(), std::invalid_argument);

    const auto explicit_full_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "mgc",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
        .resource_policy = {.cutoff_atom_threshold = 2}});
    CHECK(explicit_full_result.calculated());
    CHECK(explicit_full_result.charges->method_id() == std::string_view{"mgc"});
    CHECK(explicit_full_result.execution_policy.has_value());
    CHECK(explicit_full_result.execution_policy->mode() == calculation::ExecutionMode::full);
    CHECK(explicit_full_result.execution_issues.size() == 1);
    CHECK(explicit_full_result.execution_issues[0].kind ==
          methods::ExecutionIssueKind::resource_threshold_exceeded);

    const auto unlimited_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "mgc",
        .resource_policy = {.cutoff_atom_threshold = std::nullopt,
                            .cover_atom_threshold = std::nullopt}});
    CHECK(unlimited_result.calculated());
    CHECK(unlimited_result.charges->method_id() == std::string_view{"mgc"});
    CHECK(unlimited_result.execution_issues.empty());

    const auto throw_fn_4 = [] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "mgc",
            .execution_selection = calculation::ExecutionSelection{
                calculation::ExecutionSelectionKind::cutoff, 8.0}}));
    };
    CHECK_THROWS_AS(throw_fn_4(), std::invalid_argument);

    const auto throw_fn_5 = [] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "missing",
            .parameter_set_id = std::nullopt}));
    };
    CHECK_THROWS_AS(throw_fn_5(), std::invalid_argument);

    const auto throw_fn_6 = [] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .parameter_set_id = "missing"}));
    };
    CHECK_THROWS_AS(throw_fn_6(), std::invalid_argument);

    const auto throw_fn_7 = [] -> void {
        static_cast<void>(calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .resource_policy = {.cutoff_atom_threshold = std::nullopt,
                                .cover_atom_threshold = 10}}));
    };
    CHECK_THROWS_AS(throw_fn_7(), std::invalid_argument);
    const auto throw_fn_8 = [] -> void {
        static_cast<void>(calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .resource_policy = {.cutoff_atom_threshold = 20, .cover_atom_threshold = 10}}));
    };
    CHECK_THROWS_AS(throw_fn_8(), std::invalid_argument);

    // Explicit unsupported execution does not fall back to another mode and requires an executable
    // plan when passed to the facade.
    auto unsupported_cover_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "formal",
        .execution_selection = calculation::ExecutionSelection{
            calculation::ExecutionSelectionKind::cover, calculation::minimum_reduced_radius}});
    CHECK(!unsupported_cover_assessment.executable());
    CHECK(unsupported_cover_assessment.requires_executable_plan());
    CHECK(unsupported_cover_assessment.applicability().applicable.size() == 1);
    const auto throw_fn_9 = [&] -> void {
        static_cast<void>(calculation::calculate(std::move(unsupported_cover_assessment)));
    };
    CHECK_THROWS_AS(throw_fn_9(), std::invalid_argument);
}

TEST_CASE("calculation preserves empty-input cardinality and assessment ownership",
          "[calculation][calculation]") {

    // Empty collections and empty targets retain their distinct cardinalities: no collection
    // entries produce no assignments, while an empty molecule still produces one source assignment.
    const auto empty_collection_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector<core::Molecule>{}},
        .method_id = "formal"});
    CHECK(empty_collection_result.calculated());
    CHECK(empty_collection_result.charges->empty());

    const auto empty_target_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{core::Molecule{{}, {}, {}, "empty"}}},
        .method_id = "formal"});
    CHECK(empty_target_result.calculated());
    CHECK(empty_target_result.charges->size() == 1);
    CHECK(empty_target_result.charges->assignment(0).charges.empty());

    // AssessmentRequest rejects duplicate parameter-set IDs before filtering or applicability. This
    // remains true when the duplicates target the same or different methods and for both overloads.
    for (const auto explicit_selection : {false, true}) {
        for (const auto [first_method_id, second_method_id] :
             {std::pair{"qeq", "qeq"}, std::pair{"qeq", "eem"}}) {
            const auto lvalue_request = make_duplicate_parameter_request(
                first_method_id, second_method_id, explicit_selection);
            const auto throw_fn_10 = [&lvalue_request] -> void {
                static_cast<void>(calculation::assess(lvalue_request));
            };
            CHECK_THROWS_AS(throw_fn_10(), std::invalid_argument);

            auto rvalue_request = make_duplicate_parameter_request(
                first_method_id, second_method_id, explicit_selection);
            const auto throw_fn_11 = [&rvalue_request] -> void {
                static_cast<void>(calculation::assess(std::move(rvalue_request)));
            };
            CHECK_THROWS_AS(throw_fn_11(), std::invalid_argument);
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
    CHECK(lvalue_assessment.applicability().selected_candidate_index ==
          rvalue_assessment.applicability().selected_candidate_index);
    CHECK(lvalue_assessment.execution_policy()->mode() ==
          rvalue_assessment.execution_policy()->mode());
    const auto lvalue_result = calculation::calculate(std::move(lvalue_assessment));
    const auto rvalue_result = calculation::calculate(std::move(rvalue_assessment));
    CHECK(lvalue_result.calculated());
    CHECK(rvalue_result.calculated());
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
    CHECK(parallel_ordered_result.calculated());
    CHECK(parallel_ordered_result.charges->size() == 3);
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
