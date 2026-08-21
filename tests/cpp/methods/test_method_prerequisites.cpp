#include "support/test_molecules.h"

#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/methods/method_registry.h>

#include <cassert>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

class CoordinatesMethod final : public methods::Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{.id = "coordinates-test",
                                                          .name = "Coordinates test",
                                                          .full_name = "Coordinates test",
                                                          .publication = std::nullopt,
                                                          .priority = 0};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        auto requirements = methods::MethodRequirements{};
        requirements.coordinates = true;
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

class DenseMethod final : public methods::Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{.id = "dense-test",
                                                          .name = "Dense test",
                                                          .full_name = "Dense test",
                                                          .publication = std::nullopt,
                                                          .priority = 0};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        auto requirements = methods::MethodRequirements{};
        requirements.resources.time = methods::ComplexityTerm::atoms_cubed;
        requirements.resources.memory = methods::ComplexityTerm::atoms_squared;
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

} // namespace

auto main() -> int {
    const auto water = chargefw::test::make_water();
    const auto charged_pair = chargefw::test::make_formally_charged_pair();

    const features::PreparedMolecule prepared_water{water};
    const features::PreparedMolecule prepared_charged_pair{charged_pair};

    const auto empty_options = methods::MethodOptions{};

    const auto& registry = methods::method_registry();
    const auto* abeem = registry.find("abeem");
    const auto* dummy = registry.find("dummy");
    const auto* formal = registry.find("formal");

    assert(abeem != nullptr);
    assert(dummy != nullptr);
    assert(formal != nullptr);

    assert(dummy->check_method_prerequisites(
        {.prepared_molecule = prepared_water, .method_options = empty_options}));

    assert(formal->check_method_prerequisites(
        {.prepared_molecule = prepared_water, .method_options = empty_options}));

    const CoordinatesMethod coordinates_method;

    assert(coordinates_method.check_method_prerequisites(
        {.prepared_molecule = prepared_water, .method_options = empty_options}));

    const auto missing_coordinates = coordinates_method.check_method_prerequisites(
        {.prepared_molecule = prepared_charged_pair, .method_options = empty_options});

    assert(!missing_coordinates);
    assert(missing_coordinates.issues().size() == 1);
    assert(missing_coordinates.issues()[0].kind == methods::PrerequisiteIssueKind::missing_feature);

    const DenseMethod dense_method;

    const auto dense_result = dense_method.check_method_prerequisites(
        {.prepared_molecule = prepared_water, .method_options = empty_options});

    assert(dense_result);
    assert(dense_result.issues().empty());

    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared_collection{collection};

    const auto dummy_collection_result =
        methods::check_method_prerequisites(*dummy, prepared_collection, empty_options);

    assert(dummy_collection_result);

    const auto coordinates_collection_result =
        methods::check_method_prerequisites(coordinates_method, prepared_collection, empty_options);

    assert(!coordinates_collection_result);
    assert(coordinates_collection_result.issues().size() == 1);
    assert(coordinates_collection_result.issues()[0].kind ==
           methods::PrerequisiteIssueKind::missing_feature);

    const core::Molecule berkelium_molecule{
        std::vector{core::Atom{97}}, {}, std::vector{core::Conformer{{core::Position{}}}}};
    const features::PreparedMolecule prepared_berkelium{berkelium_molecule};

    const auto abeem_result = abeem->check_method_prerequisites(
        {.prepared_molecule = prepared_berkelium, .method_options = empty_options});
    assert(!abeem_result);
    assert(abeem_result.issues().size() == 1);
    assert(abeem_result.issues()[0].kind == methods::PrerequisiteIssueKind::unsupported_molecule);
    assert(abeem_result.issues()[0].atom_index == 0);

    return 0;
}
