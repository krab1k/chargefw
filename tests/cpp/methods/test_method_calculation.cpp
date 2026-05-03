#include "support/test_molecules.h"

#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cassert>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

class WrongSizeMethod final : public methods::Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{.id = "wrong-size-test",
                                                          .name = "Wrong size test",
                                                          .full_name = "Wrong size test",
                                                          .publication = std::nullopt,
                                                          .priority = 0};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        return {};
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput&) const
        -> charges::AtomicCharges override {
        return charges::AtomicCharges{std::vector{0.0}};
    }
};

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

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> charges::AtomicCharges override {
        const auto& molecule = input.molecule();
        const auto value = input.parameters().atom("value");

        std::vector<double> values;
        values.reserve(molecule.atom_count());

        for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
            values.push_back(value[atom_index]);
        }

        return charges::AtomicCharges{std::move(values)};
    }
};

auto atom_key(const int atomic_number,
              const parameters::AtomParameterClassificationKind classification, std::string type)
    -> parameters::AtomParameterKey {
    return {
        .atomic_number = atomic_number, .classification = classification, .type = std::move(type)};
}

auto make_collection() -> core::MoleculeCollection {
    std::vector molecules{chargefw::test::make_water(),
                          chargefw::test::make_formally_charged_pair()};

    return core::MoleculeCollection{std::move(molecules), "test-collection"};
}

auto make_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-parameters", .method_id = "atom-parameter-test", .name = "Test parameters"},
        {},
        parameters::AtomParameters{
            {{.key = atom_key(1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key =
                  atom_key(8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}},
             {.key = atom_key(7, parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "value", .value = 3.0}}},
             {.key = atom_key(17, parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "value", .value = 4.0}}}}}};
}

auto make_parameterized_candidate(const AtomParameterMethod& method,
                                  const parameters::ParameterSet& parameter_set)
    -> methods::ApplicableMethod {
    return methods::ApplicableMethod{
        .method = &method,
        .parameter_set = &parameter_set,
        .method_options = methods::MethodOptions{},
        .classifications = {
            parameters::ParameterClassification{
                parameters::AtomParameterClassification{std::vector<std::size_t>{1, 0, 0}}},
            parameters::ParameterClassification{
                parameters::AtomParameterClassification{std::vector<std::size_t>{2, 3}}}}};
}

} // namespace

auto main() -> int {
    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared{collection};

    const auto& registry = methods::method_registry();
    const auto* dummy = registry.find("dummy");

    assert(dummy != nullptr);

    const methods::ApplicableMethod dummy_candidate{
        .method = dummy,
        .parameter_set = nullptr,
        .method_options = methods::make_default_options(dummy->option_schema()),
        .classifications = {}};

    const auto dummy_charges = methods::calculate_charges(dummy_candidate, prepared);

    assert(dummy_charges.method_id() == std::string_view{"dummy"});
    assert(!dummy_charges.parameter_set_id().has_value());
    assert(dummy_charges.size() == prepared.molecule_count());

    assert(dummy_charges.assignment(0).target.molecule_index == 0);
    assert(!dummy_charges.assignment(0).target.conformer_index.has_value());
    assert(dummy_charges.assignment(0).charges.size() == collection[0].atom_count());

    for (const auto value : dummy_charges.assignment(0).charges.values()) {
        assert(value == 0.0);
    }

    const AtomParameterMethod parameterized_method;
    const auto parameter_set = make_parameter_set();

    const auto parameterized_candidate =
        make_parameterized_candidate(parameterized_method, parameter_set);

    const auto parameterized_charges =
        methods::calculate_charges(parameterized_candidate, prepared);

    assert(parameterized_charges.method_id() == std::string_view{"atom-parameter-test"});
    assert(parameterized_charges.parameter_set_id().has_value());
    assert(*parameterized_charges.parameter_set_id() == std::string_view{"test-parameters"});
    assert(parameterized_charges.size() == prepared.molecule_count());

    const auto& water_charges = parameterized_charges.assignment(0).charges;

    assert(water_charges.size() == 3);
    assert(water_charges[0] == 2.0);
    assert(water_charges[1] == 1.0);
    assert(water_charges[2] == 1.0);

    const auto& charged_pair_charges = parameterized_charges.assignment(1).charges;

    assert(charged_pair_charges.size() == 2);
    assert(charged_pair_charges[0] == 3.0);
    assert(charged_pair_charges[1] == 4.0);

    bool rejected_null_method = false;

    try {
        const methods::ApplicableMethod invalid_candidate{.method = nullptr,
                                                          .parameter_set = nullptr,
                                                          .method_options = {},
                                                          .classifications = {}};

        [[maybe_unused]] const auto invalid_charges =
            methods::calculate_charges(invalid_candidate, prepared);
    } catch (const std::invalid_argument&) {
        rejected_null_method = true;
    }

    assert(rejected_null_method);

    bool rejected_missing_classification = false;

    try {
        const methods::ApplicableMethod invalid_candidate{.method = &parameterized_method,
                                                          .parameter_set = &parameter_set,
                                                          .method_options = {},
                                                          .classifications = {}};

        [[maybe_unused]] const auto invalid_charges =
            methods::calculate_charges(invalid_candidate, prepared);
    } catch (const std::invalid_argument&) {
        rejected_missing_classification = true;
    }

    assert(rejected_missing_classification);

    bool rejected_wrong_charge_count = false;

    try {
        const WrongSizeMethod wrong_size_method;

        const methods::ApplicableMethod invalid_candidate{.method = &wrong_size_method,
                                                          .parameter_set = nullptr,
                                                          .method_options = {},
                                                          .classifications = {}};

        [[maybe_unused]] const auto invalid_charges =
            methods::calculate_charges(invalid_candidate, prepared);
    } catch (const std::invalid_argument&) {
        rejected_wrong_charge_count = true;
    }

    assert(rejected_wrong_charge_count);

    bool rejected_invalid_classification = false;

    try {
        const methods::ApplicableMethod invalid_candidate{
            .method = &parameterized_method,
            .parameter_set = &parameter_set,
            .method_options = {},
            .classifications = {
                parameters::ParameterClassification{
                    parameters::AtomParameterClassification{std::vector<std::size_t>{999, 0, 0}}},
                parameters::ParameterClassification{
                    parameters::AtomParameterClassification{std::vector<std::size_t>{2, 3}}}}};

        [[maybe_unused]] const auto invalid_charges =
            methods::calculate_charges(invalid_candidate, prepared);
    } catch (const std::invalid_argument&) {
        rejected_invalid_classification = true;
    }

    assert(rejected_invalid_classification);

    return 0;
}