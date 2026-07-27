#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cassert>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

class AtomParameterMethod final : public methods::Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{.id = "atom-parameter-test",
                                                          .name = "Atom parameter test",
                                                          .full_name = "Atom parameter test",
                                                          .publication = std::nullopt,
                                                          .priority = 0};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        auto requirements = methods::MethodRequirements{};
        requirements.bond_graph = true;
        requirements.atom_parameters = {"value"};
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& /* unused */) const
        -> chargefw::charges::AtomicCharges override {
        return chargefw::charges::AtomicCharges{std::vector<double>{}};
    }
};

auto make_collection() -> core::MoleculeCollection {
    std::vector molecules{chargefw::test::make_water(),
                          chargefw::test::make_formally_charged_pair()};

    return core::MoleculeCollection{std::move(molecules), "test-collection"};
}

auto make_wrong_method_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "wrong-method-parameters",
                                         .method_id = "other-method",
                                         .name = "Wrong method parameters"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(1, parameters::AtomParameterClassificationKind::PLAIN,
                                              "*"),
              .parameters = {{.name = "value", .value = 1.0}}}}}};
}

auto make_water_only_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "water-only-parameters",
                                         .method_id = "atom-parameter-test",
                                         .name = "Water only parameters"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}}}}};
}

auto make_collection_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "collection-parameters",
                                         .method_id = "atom-parameter-test",
                                         .name = "Collection parameters"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}},
             {.key = chargefw::test::atom_key(7, parameters::AtomParameterClassificationKind::PLAIN,
                                              "*"),
              .parameters = {{.name = "value", .value = 3.0}}},
             {.key = chargefw::test::atom_key(
                  17, parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "value", .value = 4.0}}}}}};
}

} // namespace

auto main() -> int {
    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared_collection{collection};

    const auto& registry = methods::method_registry();
    const auto* dummy = registry.find("dummy");

    assert(dummy != nullptr);

    const AtomParameterMethod atom_parameter_method;

    const std::vector<const methods::Method*> candidate_methods{dummy, &atom_parameter_method};

    const std::vector parameter_sets{make_wrong_method_parameters(), make_water_only_parameters(),
                                     make_collection_parameters()};

    const auto result =
        methods::find_applicable_methods(prepared_collection, candidate_methods, parameter_sets);

    assert(!result.empty());
    assert(result.applicable.size() == 2);

    const auto& dummy_applicable = result.applicable[0];

    assert(dummy_applicable.method == dummy);
    assert(dummy_applicable.parameter_set == nullptr);
    assert(!dummy_applicable.uses_parameters());
    assert(dummy_applicable.classifications.empty());

    const auto& parameterized_applicable = result.applicable[1];

    assert(parameterized_applicable.method == &atom_parameter_method);
    assert(parameterized_applicable.parameter_set != nullptr);
    assert(parameterized_applicable.uses_parameters());
    assert(parameterized_applicable.parameter_set->id() ==
           std::string_view{"collection-parameters"});
    assert(parameterized_applicable.classifications.size() == prepared_collection.size());

    const auto& water_classification = parameterized_applicable.classifications[0];

    assert(water_classification.atom().size() == 3);
    assert(water_classification.atom()[0] == 1);
    assert(water_classification.atom()[1] == 0);
    assert(water_classification.atom()[2] == 0);

    const auto& charged_pair_classification = parameterized_applicable.classifications[1];

    assert(charged_pair_classification.atom().size() == 2);
    assert(charged_pair_classification.atom()[0] == 2);
    assert(charged_pair_classification.atom()[1] == 3);

    assert(result.rejected.size() == 2);

    assert(result.rejected[0].method_index == 1);
    assert(result.rejected[0].parameter_set_index == std::optional<std::size_t>{0});
    assert(!result.rejected[0].issues.empty());
    assert(result.rejected[0].issues[0].kind == methods::PrerequisiteIssueKind::missing_parameters);

    assert(result.rejected[1].method_index == 1);
    assert(result.rejected[1].parameter_set_index == std::optional<std::size_t>{1});
    assert(!result.rejected[1].issues.empty());
    assert(result.rejected[1].issues[0].kind ==
           methods::PrerequisiteIssueKind::parameter_classification_failed);

    return 0;
}
