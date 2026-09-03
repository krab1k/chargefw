#include "support/test_molecules.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/methods/method_registry.h>

#include <limits>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <utility>
#include <vector>

namespace methods = chargefw::methods;

TEST_CASE("prerequisite issue kinds convert to stable strings", "[methods][prerequisites]") {
    CHECK(methods::to_string(methods::PrerequisiteIssueKind::invalid_options) == "invalid_options");
    CHECK(methods::to_string(methods::PrerequisiteIssueKind::missing_parameters) ==
          "missing_parameters");
}

namespace core = chargefw::core;
namespace features = chargefw::features;

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

TEST_CASE("method prerequisites detect missing features and unsupported molecules",
          "[methods][method-prerequisites]") {
    const auto water = chargefw::test::make_water();
    const auto charged_pair = chargefw::test::make_formally_charged_pair();

    const features::PreparedMolecule prepared_water{water};
    const features::PreparedMolecule prepared_charged_pair{charged_pair};

    const auto empty_options = methods::MethodOptions{};

    const auto& registry = methods::method_registry();
    const auto* abeem = registry.find("abeem");
    const auto* dummy = registry.find("dummy");
    const auto* formal = registry.find("formal");

    REQUIRE(abeem != nullptr);
    REQUIRE(dummy != nullptr);
    REQUIRE(formal != nullptr);

    CHECK(dummy->check_method_prerequisites(
        {.prepared_molecule = prepared_water, .method_options = empty_options}));

    CHECK(formal->check_method_prerequisites(
        {.prepared_molecule = prepared_water, .method_options = empty_options}));

    const CoordinatesMethod coordinates_method;

    CHECK(coordinates_method.check_method_prerequisites(
        {.prepared_molecule = prepared_water, .method_options = empty_options}));

    const auto missing_coordinates = coordinates_method.check_method_prerequisites(
        {.prepared_molecule = prepared_charged_pair, .method_options = empty_options});

    CHECK(!missing_coordinates);
    REQUIRE(missing_coordinates.issues().size() == 1);
    CHECK(missing_coordinates.issues()[0].kind == methods::PrerequisiteIssueKind::missing_feature);

    const core::Molecule coincident_atoms{
        std::vector{core::Atom{1}, core::Atom{1}},
        {},
        std::vector{core::Conformer{{core::Position{.x = 1.5, .y = -2.25, .z = 3.75},
                                     core::Position{.x = 1.5, .y = -2.25, .z = 3.75}},
                                    "first"},
                    core::Conformer{{core::Position{}, core::Position{.x = -0.0}}, "second"}}};
    const features::PreparedMolecule prepared_coincident_atoms{coincident_atoms};
    const auto coincident_atoms_result = coordinates_method.check_method_prerequisites(
        {.prepared_molecule = prepared_coincident_atoms, .method_options = empty_options});
    CHECK(!coincident_atoms_result);
    REQUIRE(coincident_atoms_result.issues().size() == 2);
    CHECK(coincident_atoms_result.issues()[0].kind ==
          methods::PrerequisiteIssueKind::invalid_geometry);
    CHECK(coincident_atoms_result.issues()[0].atom_index == 1);
    CHECK(coincident_atoms_result.issues()[0].conformer_index == 0);
    CHECK(coincident_atoms_result.issues()[0].message.contains("conformer 1"));
    CHECK(coincident_atoms_result.issues()[0].message.contains("atom 1 (H, formal charge 0)"));
    CHECK(coincident_atoms_result.issues()[1].kind ==
          methods::PrerequisiteIssueKind::invalid_geometry);
    CHECK(coincident_atoms_result.issues()[1].atom_index == 1);
    CHECK(coincident_atoms_result.issues()[1].conformer_index == 1);
    CHECK(coincident_atoms_result.issues()[1].message.contains("conformer 2"));

    const core::Molecule nonfinite_atom{
        std::vector{core::Atom{1}},
        {},
        std::vector{core::Conformer{{core::Position{.x = std::numeric_limits<double>::quiet_NaN()}},
                                    "nonfinite"}}};
    const features::PreparedMolecule prepared_nonfinite_atom{nonfinite_atom};
    const auto nonfinite_atom_result = coordinates_method.check_method_prerequisites(
        {.prepared_molecule = prepared_nonfinite_atom, .method_options = empty_options});
    CHECK(!nonfinite_atom_result);
    REQUIRE(nonfinite_atom_result.issues().size() == 1);
    CHECK(nonfinite_atom_result.issues()[0].kind ==
          methods::PrerequisiteIssueKind::invalid_geometry);
    CHECK(nonfinite_atom_result.issues()[0].atom_index == 0);
    CHECK(nonfinite_atom_result.issues()[0].conformer_index == 0);
    CHECK(nonfinite_atom_result.issues()[0].message ==
          "method 'coordinates-test', conformer 1: atom 1 (H, formal charge 0) has non-finite "
          "coordinates");

    const DenseMethod dense_method;

    const auto dense_result = dense_method.check_method_prerequisites(
        {.prepared_molecule = prepared_water, .method_options = empty_options});

    CHECK(dense_result);
    CHECK(dense_result.issues().empty());

    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared_collection{collection};

    const auto dummy_collection_result =
        methods::check_method_prerequisites(*dummy, prepared_collection, empty_options);

    CHECK(dummy_collection_result);

    const auto coordinates_collection_result =
        methods::check_method_prerequisites(coordinates_method, prepared_collection, empty_options);

    CHECK(!coordinates_collection_result);
    REQUIRE(coordinates_collection_result.issues().size() == 1);
    CHECK(coordinates_collection_result.issues()[0].kind ==
          methods::PrerequisiteIssueKind::missing_feature);
    CHECK(coordinates_collection_result.issues()[0].molecule_index == 1);
    CHECK(coordinates_collection_result.issues()[0].message.contains("molecule 2"));

    const core::Molecule berkelium_molecule{
        std::vector{core::Atom{97}}, {}, std::vector{core::Conformer{{core::Position{}}}}};
    const features::PreparedMolecule prepared_berkelium{berkelium_molecule};

    const auto abeem_result = abeem->check_method_prerequisites(
        {.prepared_molecule = prepared_berkelium, .method_options = empty_options});
    CHECK(!abeem_result);
    REQUIRE(abeem_result.issues().size() == 1);
    CHECK(abeem_result.issues()[0].kind == methods::PrerequisiteIssueKind::unsupported_molecule);
    CHECK(abeem_result.issues()[0].atom_index == 0);
    CHECK(abeem_result.issues()[0].message.contains("atom 1 (Bk, formal charge 0)"));
}
