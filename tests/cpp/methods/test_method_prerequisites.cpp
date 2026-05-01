#include "support/test_molecules.h"

#include <chargefw/methods/method.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/methods/method_options.h>

#include <cassert>
#include <span>

namespace methods = chargefw::methods;

namespace {

class CoordinatesMethod final : public methods::Method {
public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{
            .id = "coordinates-test",
            .name = "Coordinates test",
            .full_name = "Coordinates test",
            .publication = std::nullopt,
            .priority = 0
        };

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        return methods::MethodRequirements{.coordinates = true};
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput&) const
        -> chargefw::charges::AtomicCharges override {
        return chargefw::charges::AtomicCharges{{}};
    }
};

class DenseMethod final : public methods::Method {
public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{
            .id = "dense-test",
            .name = "Dense test",
            .full_name = "Dense test",
            .publication = std::nullopt,
            .priority = 0
        };

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        return methods::MethodRequirements{
            .resources = methods::ResourceRequirements{
                .time = methods::ComplexityTerm::atoms_cubed,
                .memory = methods::ComplexityTerm::atoms_squared,
                .large_molecule_atom_threshold = 2,
                .reject_large_without_reduction = true
            }
        };
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput&) const
        -> chargefw::charges::AtomicCharges override {
        return chargefw::charges::AtomicCharges{{}};
    }
};

} // namespace

auto main() -> int {
    const auto water = chargefw::test::make_water();
    const auto charged_pair = chargefw::test::make_formally_charged_pair();

    const auto empty_options = methods::MethodOptions{};

    const auto& registry = methods::method_registry();
    const auto* dummy = registry.find("dummy");
    const auto* formal = registry.find("formal");

    assert(dummy != nullptr);
    assert(formal != nullptr);

    assert(dummy->check_method_prerequisites({
        .molecule = water,
        .method_options = empty_options
    }));

    assert(formal->check_method_prerequisites({
        .molecule = water,
        .method_options = empty_options
    }));

    const CoordinatesMethod coordinates_method;

    assert(coordinates_method.check_method_prerequisites({
        .molecule = water,
        .method_options = empty_options
    }));

    const auto missing_coordinates = coordinates_method.check_method_prerequisites({
        .molecule = charged_pair,
        .method_options = empty_options
    });

    assert(!missing_coordinates);
    assert(missing_coordinates.issues().size() == 1);
    assert(missing_coordinates.issues()[0].kind ==
           methods::PrerequisiteIssueKind::missing_feature);

    const DenseMethod dense_method;

    const auto dense_result = dense_method.check_method_prerequisites({
        .molecule = water,
        .method_options = empty_options
    });

    assert(!dense_result);
    assert(dense_result.issues().size() == 1);
    assert(dense_result.issues()[0].kind ==
           methods::PrerequisiteIssueKind::resource_limit);

    return 0;
}