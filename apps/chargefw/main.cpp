#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/features/topology_features.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/methods/parameter_prerequisites.h>
#include <chargefw/parameters/atom_parameters.h>
#include <chargefw/parameters/parameter_key.h>
#include <chargefw/parameters/parameter_set.h>
#include <chargefw/parameters/parameter_set_metadata.h>

#include <iostream>
#include <memory>
#include <optional>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

auto make_water() -> core::Molecule {
    std::vector atoms{core::Atom{8, 0, "O"}, core::Atom{1, 0, "H1"}, core::Atom{1, 0, "H2"}};

    std::vector bonds{core::Bond{0, 1, core::BondOrder::SINGLE},
                      core::Bond{0, 2, core::BondOrder::SINGLE}};

    std::vector positions{core::Position{.x = 0.0000, .y = 0.0000, .z = 0.0000},
                          core::Position{.x = 0.9572, .y = 0.0000, .z = 0.0000},
                          core::Position{.x = -0.2390, .y = 0.9270, .z = 0.0000}};

    std::vector conformers{core::Conformer{std::move(positions), "model-1"}};

    return core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers), "water"};
}

auto atom_key(const int atomic_number,
              const parameters::AtomParameterClassificationKind classification, std::string type)
    -> parameters::AtomParameterKey {
    return {
        .atomic_number = atomic_number, .classification = classification, .type = std::move(type)};
}

auto make_matching_water_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "demo-water-parameters",
                                         .method_id = "demo-atom-parameter-method",
                                         .name = "Demo water parameters"},
        {},
        parameters::AtomParameters{
            {{.key = atom_key(1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key =
                  atom_key(8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}}}}};
}

auto make_nonmatching_water_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "demo-incomplete-parameters",
                                         .method_id = "demo-atom-parameter-method",
                                         .name = "Incomplete demo parameters"},
        {},
        parameters::AtomParameters{
            {{.key = atom_key(1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}}}}};
}

class DemoAtomParameterMethod final : public methods::Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        static constexpr methods::MethodMetadata metadata{
            .id = "demo-atom-parameter-method",
            .name = "Demo atom parameter method",
            .full_name = "Demo method requiring atom parameters",
            .publication = std::nullopt,
            .priority = 0};

        return metadata;
    }

    [[nodiscard]] auto requirements() const noexcept -> methods::MethodRequirements override {
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
        return charges::AtomicCharges{std::vector<double>(input.molecule.atom_count(), 0.0)};
    }
};

auto print_result(std::string_view label, const methods::PrerequisiteResult& result) -> void {
    std::cout << label << ": " << (result ? "ok" : "not ok") << '\n';

    for (const auto& issue : result.issues()) {
        std::cout << "  - " << issue.message << '\n';
    }
}

auto print_parameter_result(std::string_view label,
                            const methods::ParameterPrerequisiteResult& result) -> void {
    std::cout << label << ": " << (result ? "ok" : "not ok") << '\n';

    for (const auto& issue : result.issues) {
        std::cout << "  - " << issue.message << '\n';
    }

    if (result.classification.has_value()) {
        std::cout << "  atom classification:";
        for (const auto index : result.classification->atom().parameter_entry_indices()) {
            std::cout << ' ' << index;
        }
        std::cout << '\n';
    }
}

auto calculate_and_print(const methods::Method& method, const core::Molecule& molecule,
                         const features::TopologyFeatures& topology) -> void {
    const auto method_options = methods::make_default_options(method.option_schema());

    const methods::CalculationInput input{.molecule = molecule,
                                          .topology = topology,
                                          .geometry = nullptr,
                                          .method_options = method_options};

    const auto charges = method.calculate(input);

    std::cout << method.id() << " charges:";
    for (const auto charge : charges.values()) {
        std::cout << ' ' << charge;
    }
    std::cout << "  total=" << charges.total() << '\n';
}

} // namespace

auto main() -> int {
    const auto water = make_water();
    const features::TopologyFeatures topology{water};

    const auto& registry = methods::method_registry();

    const auto* dummy = registry.find("dummy");
    const auto* formal = registry.find("formal");

    if (dummy == nullptr || formal == nullptr) {
        std::cerr << "Required built-in methods are missing.\n";
        return 1;
    }

    const auto dummy_options = methods::make_default_options(dummy->option_schema());
    const auto formal_options = methods::make_default_options(formal->option_schema());

    print_result(
        "dummy method prerequisites",
        dummy->check_method_prerequisites({.molecule = water, .method_options = dummy_options}));

    print_result(
        "formal method prerequisites",
        formal->check_method_prerequisites({.molecule = water, .method_options = formal_options}));

    calculate_and_print(*dummy, water, topology);
    calculate_and_print(*formal, water, topology);

    const DemoAtomParameterMethod demo_method;
    const auto demo_options = methods::make_default_options(demo_method.option_schema());

    print_result("demo parameterized method prerequisites",
                 demo_method.check_method_prerequisites(
                     {.molecule = water, .method_options = demo_options}));

    const auto matching_parameters = make_matching_water_parameters();
    const auto nonmatching_parameters = make_nonmatching_water_parameters();

    print_parameter_result(
        "matching parameter prerequisites",
        methods::check_parameter_prerequisites(demo_method, {.molecule = water,
                                                             .topology = topology,
                                                             .parameter_set = matching_parameters,
                                                             .classification_options = {}}));

    print_parameter_result("nonmatching parameter prerequisites",
                           methods::check_parameter_prerequisites(
                               demo_method, {.molecule = water,
                                             .topology = topology,
                                             .parameter_set = nonmatching_parameters,
                                             .classification_options = {}}));

    return 0;
}