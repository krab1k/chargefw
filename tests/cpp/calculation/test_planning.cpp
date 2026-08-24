#include "support/test_molecules.h"
#include "support/test_parameters.h"

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

#include <algorithm>
#include <cstdint>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

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

TEST_CASE("planning selects the highest-priority full execution plan", "[calculation][planning]") {
    const auto prepared = make_prepared_water();
    const FixedChargeMethod higher_priority{"higher", 10, 10.0};
    const FixedChargeMethod lower_priority{"lower", 1, 1.0};
    const std::vector<const methods::Method*> methods{&lower_priority, &higher_priority};
    const std::vector<chargefw::parameters::ParameterSet> parameters;

    const auto applicability = methods::find_applicable_methods(
        {.molecules = prepared, .methods = methods, .parameter_sets = parameters});
    const auto* selected = calculation::select_applicable_method(applicability);
    REQUIRE(selected != nullptr);
    const auto plan =
        calculation::select_execution_plan(applicability, calculation::ExecutionSelection{});
    REQUIRE(plan.has_value());
    CHECK(plan->selected == selected);
    CHECK(plan->policy.mode() == calculation::ExecutionMode::full);
    CHECK(plan->issues.empty());
    CHECK(applicability.applicable.size() == 2);
}

TEST_CASE("planning ranks automatic candidates and supports permissive classification",
          "[calculation][planning]") {
    const auto prepared = make_prepared_water();
    const std::vector<chargefw::parameters::ParameterSet> parameters;

    const FixedChargeMethod alpha{"alpha", 1, 2.0};
    const FixedChargeMethod beta{"beta", 1, 3.0};
    const std::vector<const methods::Method*> tied_methods{&beta, &alpha};
    const auto tied_result = calculate_automatically(prepared, tied_methods, parameters);
    REQUIRE(tied_result.charges.has_value());
    CHECK(tied_result.charges->method_id() == std::string_view{"alpha"});

    const ParameterizedFixedChargeMethod parameterized{"parameterized", 0, 4.0};
    const std::vector<const methods::Method*> parameterized_methods{&parameterized};
    const std::vector parameter_sets{make_parameter_set("alpha", "parameterized", 1),
                                     make_parameter_set("zeta", "parameterized", 10)};
    const auto parameterized_result =
        calculate_automatically(prepared, parameterized_methods, parameter_sets);
    REQUIRE(parameterized_result.charges.has_value());
    CHECK(parameterized_result.charges->parameter_set_id() == std::string_view{"zeta"});

    const std::vector<const methods::Method*> no_methods;
    const auto no_result = calculate_automatically(prepared, no_methods, parameters);
    CHECK_FALSE(no_result.calculated());
    CHECK(no_result.applicability.empty());

    const core::MoleculeCollection double_bonded_collection{
        std::vector{make_double_bonded_carbons()}};
    const features::PreparedMoleculeCollection double_bonded_prepared{double_bonded_collection};
    const ParameterValueMethod permissive_method;
    const std::vector<const methods::Method*> permissive_methods{&permissive_method};
    const std::vector permissive_parameter_sets{make_permissive_parameter_set()};

    const auto strict_result = calculate_automatically(double_bonded_prepared, permissive_methods,
                                                       permissive_parameter_sets);
    CHECK_FALSE(strict_result.calculated());

    const auto permissive_result =
        calculate_automatically(double_bonded_prepared, permissive_methods,
                                permissive_parameter_sets, {.permissive_types = true});
    REQUIRE(permissive_result.charges.has_value());
    CHECK(permissive_result.applicability.applicable.size() == 1);
    CHECK(permissive_result.charges->assignment(0).charges[0] == 3.0);
    CHECK(permissive_result.charges->assignment(0).charges[1] == 3.0);
}

TEST_CASE("explicit unsupported execution has no selected plan or fallback",
          "[calculation][planning]") {
    auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "formal",
        .execution_selection = calculation::ExecutionSelection{
            calculation::ExecutionSelectionKind::cover, calculation::minimum_reduced_radius}});

    CHECK_FALSE(assessment.executable());
    CHECK(assessment.requires_executable_plan());
    CHECK_FALSE(assessment.applicability().selected_candidate_index.has_value());
    CHECK_FALSE(assessment.execution_policy().has_value());
    CHECK(assessment.execution_issues().empty());
    REQUIRE(assessment.applicability().applicable.size() == 1);

    const auto& cover_assessment = assessment.applicability().applicable[0].execution_assessments;
    const auto cover = std::ranges::find_if(cover_assessment, [](const auto& value) {
        return value.mode == calculation::ExecutionMode::cover;
    });
    REQUIRE(cover != cover_assessment.end());
    CHECK(cover->availability == methods::ExecutionAvailability::unsupported);
    REQUIRE(cover->issues.size() == 1);
    CHECK(cover->issues[0].kind == methods::ExecutionIssueKind::unsupported_execution_mode);

    const auto calculate_unsupported_cover = [&] {
        static_cast<void>(calculation::calculate(std::move(assessment)));
    };
    CHECK_THROWS_AS(calculate_unsupported_cover(), std::invalid_argument);
}
