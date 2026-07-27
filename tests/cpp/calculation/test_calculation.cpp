#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <cassert>
#include <chargefw/calculation/calculation.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <optional>
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

template <typename Callable> auto throws_invalid_argument(Callable&& callable) -> bool {
    try {
        std::forward<Callable>(callable)();
    } catch (const std::invalid_argument&) {
        return true;
    }

    return false;
}

} // namespace

auto main() -> int {
    const auto prepared = make_prepared_water();
    const FixedChargeMethod higher_priority{"higher", 10, 10.0};
    const FixedChargeMethod lower_priority{"lower", 1, 1.0};
    const std::vector<const methods::Method*> methods{&lower_priority, &higher_priority};
    const std::vector<chargefw::parameters::ParameterSet> parameters;

    const auto result = calculation::calculate(calculation::CalculationRequest{
        .molecules = prepared, .candidate_methods = methods, .parameter_sets = parameters});

    assert(result.calculated());
    assert(result.charges->method_id() == std::string_view{"higher"});
    assert(result.charges->size() == 1);
    assert(result.charges->assignment(0).charges[0] == 10.0);
    assert(result.applicability.applicable.size() == 2);

    const FixedChargeMethod alpha{"alpha", 1, 2.0};
    const FixedChargeMethod beta{"beta", 1, 3.0};
    const std::vector<const methods::Method*> tied_methods{&beta, &alpha};

    const auto tied_result = calculation::calculate(calculation::CalculationRequest{
        .molecules = prepared, .candidate_methods = tied_methods, .parameter_sets = parameters});

    assert(tied_result.calculated());
    assert(tied_result.charges->method_id() == std::string_view{"alpha"});
    assert(tied_result.charges->assignment(0).charges[0] == 2.0);

    const ParameterizedFixedChargeMethod parameterized{"parameterized", 0, 4.0};
    const std::vector<const methods::Method*> parameterized_methods{&parameterized};
    const std::vector parameter_sets{make_parameter_set("alpha", "parameterized", 1),
                                     make_parameter_set("zeta", "parameterized", 10)};

    const auto parameterized_result = calculation::calculate(
        calculation::CalculationRequest{.molecules = prepared,
                                        .candidate_methods = parameterized_methods,
                                        .parameter_sets = parameter_sets});

    assert(parameterized_result.calculated());
    assert(parameterized_result.charges->method_id() == std::string_view{"parameterized"});
    assert(parameterized_result.charges->parameter_set_id() == std::string_view{"zeta"});

    const std::vector<const methods::Method*> no_methods;
    const auto no_result = calculation::calculate(calculation::CalculationRequest{
        .molecules = prepared, .candidate_methods = no_methods, .parameter_sets = parameters});

    assert(!no_result.calculated());
    assert(no_result.applicability.empty());

    const auto application_result =
        calculation::calculate(calculation::ApplicationCalculationRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .parameter_set_id = std::nullopt});

    assert(application_result.calculated());
    assert(application_result.charges->method_id() == std::string_view{"formal"});
    assert(application_result.charges->size() == 1);

    assert(throws_invalid_argument([] {
        static_cast<void>(calculation::calculate(calculation::ApplicationCalculationRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "missing",
            .parameter_set_id = std::nullopt}));
    }));

    assert(throws_invalid_argument([] {
        static_cast<void>(calculation::calculate(calculation::ApplicationCalculationRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {},
            .method_id = "formal",
            .parameter_set_id = "missing"}));
    }));

    return 0;
}
