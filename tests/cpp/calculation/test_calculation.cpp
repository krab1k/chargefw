#include "support/test_assertions.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <algorithm>
#include <cassert>
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

auto main() -> int {
    static_assert(!HasPublicParameterSets<calculation::AssessmentResult>);
    static_assert(!HasPublicSelectedCandidate<calculation::AssessmentResult>);
    static_assert(std::is_move_constructible_v<calculation::AssessmentResult>);
    static_assert(!std::is_move_assignable_v<calculation::AssessmentResult>);

    const auto prepared = make_prepared_water();
    const FixedChargeMethod higher_priority{"higher", 10, 10.0};
    const FixedChargeMethod lower_priority{"lower", 1, 1.0};
    const std::vector<const methods::Method*> methods{&lower_priority, &higher_priority};
    const std::vector<chargefw::parameters::ParameterSet> parameters;

    const auto applicability = methods::find_applicable_methods(
        {.molecules = prepared, .methods = methods, .parameter_sets = parameters});
    const auto* selected = calculation::select_applicable_method(applicability);
    assert(selected != nullptr);
    const auto plan =
        calculation::select_execution_plan(applicability, calculation::ExecutionSelection{});
    assert(plan.has_value());
    assert(plan->selected == selected);
    assert(plan->policy.mode() == calculation::ExecutionMode::full);
    assert(plan->issues.empty());
    const auto result = calculation::calculate({.molecules = prepared, .selected = *selected});
    assert(result.charges.method_id() == std::string_view{"higher"});
    assert(result.charges.size() == 1);
    assert(result.charges.assignment(0).charges[0] == 10.0);
    assert(applicability.applicable.size() == 2);

    const auto parallel_collection = core::MoleculeCollection{
        std::vector{chargefw::test::make_water(), chargefw::test::make_water()}};
    const features::PreparedMoleculeCollection parallel_prepared{parallel_collection};
    const auto parallel_result =
        calculation::calculate({.molecules = parallel_prepared,
                                .selected = *selected,
                                .execution_policy = calculation::ExecutionPolicy{},
                                .max_threads = 2});
    assert(parallel_result.charges.size() == 2);
    assert(parallel_result.charges.assignment(0).target.molecule_index == 0);
    assert(parallel_result.charges.assignment(1).target.molecule_index == 1);
    assert(parallel_result.charges.assignment(0).charges[0] == 10.0);
    assert(parallel_result.charges.assignment(1).charges[0] == 10.0);

    assert(chargefw::test::throws_invalid_argument([&] -> void {
        static_cast<void>(
            calculation::calculate({.molecules = prepared,
                                    .selected = *selected,
                                    .execution_policy = calculation::ExecutionPolicy{},
                                    .max_threads = std::numeric_limits<std::size_t>::max()}));
    }));

    const FixedChargeMethod alpha{"alpha", 1, 2.0};
    const FixedChargeMethod beta{"beta", 1, 3.0};
    const std::vector<const methods::Method*> tied_methods{&beta, &alpha};

    const auto tied_result = calculate_automatically(prepared, tied_methods, parameters);

    if (!tied_result.charges.has_value()) {
        return 1;
    }
    assert(tied_result.charges->method_id() == std::string_view{"alpha"});
    assert(tied_result.charges->assignment(0).charges[0] == 2.0);

    const ParameterizedFixedChargeMethod parameterized{"parameterized", 0, 4.0};
    const std::vector<const methods::Method*> parameterized_methods{&parameterized};
    const std::vector parameter_sets{make_parameter_set("alpha", "parameterized", 1),
                                     make_parameter_set("zeta", "parameterized", 10)};

    const auto parameterized_result =
        calculate_automatically(prepared, parameterized_methods, parameter_sets);

    if (!parameterized_result.charges.has_value()) {
        return 1;
    }
    assert(parameterized_result.charges->method_id() == std::string_view{"parameterized"});
    assert(parameterized_result.charges->parameter_set_id() == std::string_view{"zeta"});

    const std::vector<const methods::Method*> no_methods;
    const auto no_result = calculate_automatically(prepared, no_methods, parameters);

    assert(!no_result.calculated());
    assert(no_result.applicability.empty());

    const core::MoleculeCollection double_bonded_collection{
        std::vector{make_double_bonded_carbons()}};
    const features::PreparedMoleculeCollection double_bonded_prepared{double_bonded_collection};
    const ParameterValueMethod permissive_method;
    const std::vector<const methods::Method*> permissive_methods{&permissive_method};
    const std::vector permissive_parameter_sets{make_permissive_parameter_set()};

    const auto strict_permissive_result = calculate_automatically(
        double_bonded_prepared, permissive_methods, permissive_parameter_sets);
    assert(!strict_permissive_result.calculated());

    const auto permissive_calculation_result =
        calculate_automatically(double_bonded_prepared, permissive_methods,
                                permissive_parameter_sets, {.permissive_types = true});
    assert(permissive_calculation_result.calculated());
    assert(permissive_calculation_result.applicability.applicable.size() == 1);
    assert(permissive_calculation_result.charges->assignment(0).charges[0] == 3.0);
    assert(permissive_calculation_result.charges->assignment(0).charges[1] == 3.0);

    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
            .parameter_sets = {make_permissive_peoe_parameter_set()},
            .method_id = "peoe",
            .parameter_set_id = "permissive-peoe-parameters"}));
    }));

    const auto permissive_application_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
        .parameter_sets = {make_permissive_peoe_parameter_set()},
        .method_id = "peoe",
        .parameter_set_id = "permissive-peoe-parameters",
        .classification_options = {.permissive_types = true}});
    assert(permissive_application_result.calculated());
    assert(permissive_application_result.charges->method_id() == std::string_view{"peoe"});

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
    assert(owned_parameterized_result.calculated());
    assert(owned_parameterized_result.applicability.applicable.size() == 1);
    assert(owned_parameterized_result.applicability.applicable[0].method_id == "peoe");
    assert(owned_parameterized_result.applicability.applicable[0].parameter_set_id ==
           std::optional<std::string>{"permissive-peoe-parameters"});
    assert(owned_parameterized_result.applicability.selected_candidate_index ==
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
    assert(consuming_assessment.executable());
    assert(consuming_assessment.applicability().selected_candidate_index ==
           std::optional<std::size_t>{0});
    const auto consuming_result = calculation::calculate(std::move(consuming_assessment));
    assert(consuming_result.calculated());
    assert(consuming_result.charges->size() == 2);
    assert(consuming_result.charges->assignment(0).target.molecule_index == 0);
    assert(consuming_result.charges->assignment(1).target.molecule_index == 1);

    auto peoe_options = methods::MethodOptions{};
    peoe_options.set("iters", 1);
    const auto configured_peoe_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
        .parameter_sets = {make_permissive_peoe_parameter_set()},
        .method_id = "peoe",
        .parameter_set_id = "permissive-peoe-parameters",
        .method_options = {{"peoe", peoe_options}},
        .classification_options = {.permissive_types = true}});
    assert(configured_peoe_result.calculated());
    assert(configured_peoe_result.effective_method_options.has_value());
    assert(configured_peoe_result.effective_method_options->get<int>("iters") == 1);

    auto invalid_peoe_options = methods::MethodOptions{};
    invalid_peoe_options.set("iters", std::string{"one"});
    assert(chargefw::test::throws_invalid_argument([&] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
            .parameter_sets = {make_permissive_peoe_parameter_set()},
            .method_id = "peoe",
            .parameter_set_id = "permissive-peoe-parameters",
            .method_options = {{"peoe", invalid_peoe_options}},
            .classification_options = {.permissive_types = true}}));
    }));

    const auto rejected_explicit_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "smpqeq"});
    assert(!rejected_explicit_assessment.executable());
    assert(rejected_explicit_assessment.applicability().rejected.size() == 1);
    assert(rejected_explicit_assessment.applicability().rejected[0].method_id == "smpqeq");
    assert(!rejected_explicit_assessment.applicability().rejected[0].issues.empty());

    const auto rejected_parameter_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
        .parameter_sets = {make_permissive_peoe_parameter_set()},
        .method_id = "peoe",
        .parameter_set_id = "permissive-peoe-parameters"});
    assert(!rejected_parameter_assessment.executable());
    assert(rejected_parameter_assessment.applicability().rejected.size() == 1);
    assert(rejected_parameter_assessment.applicability().rejected[0].method_id == "peoe");
    assert(rejected_parameter_assessment.applicability().rejected[0].parameter_set_id ==
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
    assert(rejected_smpqeq != automatic_rejected_result.applicability.rejected.end());
    assert(!rejected_smpqeq->issues.empty());

    const auto application_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "formal",
        .parameter_set_id = std::nullopt,
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});

    if (!application_result.charges.has_value()) {
        return 1;
    }
    assert(application_result.charges->method_id() == std::string_view{"formal"});
    assert(application_result.charges->size() == 1);
    assert(application_result.execution_policy.has_value());
    assert(application_result.execution_policy->mode() == calculation::ExecutionMode::full);
    assert(application_result.execution_issues.empty());

    auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "formal",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});
    assert(assessment.executable());
    assert(assessment.applicability().applicable.size() == 1);
    assert(assessment.applicability().selected_candidate_index == std::optional<std::size_t>{0});
    assert(assessment.applicability().applicable[0].method_id == "formal");
    assert(assessment.execution_policy()->mode() == calculation::ExecutionMode::full);

    const auto assessed_result = calculation::calculate(std::move(assessment), 1);
    assert(assessed_result.calculated());
    assert(assessed_result.charges->method_id() == std::string_view{"formal"});
    assert(assessed_result.metrics.applicability_seconds >= 0.0);

    const auto automatic_fallback_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .resource_policy = {.full_atom_threshold = 2}});
    assert(automatic_fallback_result.calculated());
    assert(automatic_fallback_result.charges->method_id() == std::string_view{"eqeq"});
    assert(automatic_fallback_result.execution_policy.has_value());
    assert(automatic_fallback_result.execution_policy->mode() ==
           calculation::ExecutionMode::cutoff);
    assert(automatic_fallback_result.execution_policy->radius() ==
           std::optional<double>{calculation::default_automatic_reduced_radius});
    assert(automatic_fallback_result.execution_issues.empty());

    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "mgc",
            .resource_policy = {.full_atom_threshold = 2}}));
    }));

    const auto explicit_full_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "mgc",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
        .resource_policy = {.full_atom_threshold = 2}});
    assert(explicit_full_result.calculated());
    assert(explicit_full_result.charges->method_id() == std::string_view{"mgc"});
    assert(explicit_full_result.execution_policy.has_value());
    assert(explicit_full_result.execution_policy->mode() == calculation::ExecutionMode::full);
    assert(explicit_full_result.execution_issues.size() == 1);
    assert(explicit_full_result.execution_issues[0].kind ==
           methods::ExecutionIssueKind::resource_threshold_exceeded);

    const auto unlimited_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "mgc",
        .resource_policy = {.full_atom_threshold = std::nullopt}});
    assert(unlimited_result.calculated());
    assert(unlimited_result.charges->method_id() == std::string_view{"mgc"});
    assert(unlimited_result.execution_issues.empty());

    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "mgc",
            .execution_selection = calculation::ExecutionSelection{
                calculation::ExecutionSelectionKind::cutoff, 8.0}}));
    }));

    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "missing",
            .parameter_set_id = std::nullopt}));
    }));

    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .parameter_set_id = "missing"}));
    }));

    return 0;
}
