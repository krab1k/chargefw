#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/prepared_molecule_collection.h>
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

auto make_formally_charged_pair() -> core::Molecule {
    std::vector atoms{core::Atom{7, 1, "N"}, core::Atom{17, -1, "Cl"}};

    return core::Molecule{std::move(atoms), {}, {}, "charged-pair"};
}

auto make_demo_collection() -> core::MoleculeCollection {
    std::vector<core::Molecule> molecules;
    molecules.push_back(make_water());
    molecules.push_back(make_formally_charged_pair());

    return core::MoleculeCollection{std::move(molecules), "demo collection"};
}

auto atom_key(const int atomic_number,
              const parameters::AtomParameterClassificationKind classification, std::string type)
    -> parameters::AtomParameterKey {
    return {
        .atomic_number = atomic_number, .classification = classification, .type = std::move(type)};
}

auto make_water_parameters() -> parameters::ParameterSet {
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
        return charges::AtomicCharges{std::vector<double>(molecule.atom_count(), 0.0)};
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

auto print_collection_parameter_result(std::string_view label,
                                       const methods::CollectionParameterPrerequisiteResult& result)
    -> void {
    std::cout << label << ": " << (result ? "ok" : "not ok") << '\n';

    for (const auto& issue : result.issues) {
        std::cout << "  - " << issue.message << '\n';
    }

    if (result) {
        std::cout << "  classifications: " << result.classifications.size() << '\n';
    }
}

auto calculate_and_print(const methods::Method& method,
                         const features::PreparedMolecule& prepared_molecule) -> void {
    const auto method_options = methods::make_default_options(method.option_schema());

    const methods::CalculationInput input{prepared_molecule, method_options};

    const auto charges = method.calculate(input);

    std::cout << method.id() << " charges:";

    for (const auto charge : charges.values()) {
        std::cout << ' ' << charge;
    }

    std::cout << "  total=" << charges.total() << '\n';
}

} // namespace

auto main() -> int {
    const auto collection = make_demo_collection();
    const auto prepared = features::PreparedMoleculeCollection{collection};

    const auto& registry = methods::method_registry();

    const auto* dummy = registry.find("dummy");
    const auto* formal = registry.find("formal");
    const auto* veem = registry.find("veem");

    if (dummy == nullptr || formal == nullptr || veem == nullptr) {
        std::cerr << "Required built-in methods are missing.\n";
        return 1;
    }

    const auto dummy_options = methods::make_default_options(dummy->option_schema());
    const auto formal_options = methods::make_default_options(formal->option_schema());
    const auto veem_options = methods::make_default_options(veem->option_schema());

    print_result("dummy collection method prerequisites",
                 methods::check_method_prerequisites(*dummy, prepared, dummy_options));

    print_result("formal collection method prerequisites",
                 methods::check_method_prerequisites(*formal, prepared, formal_options));

    print_result("veem collection method prerequisites",
                 methods::check_method_prerequisites(*veem, prepared, veem_options));

    std::cout << "\nCharges for first molecule:\n";
    calculate_and_print(*dummy, prepared[0]);
    calculate_and_print(*formal, prepared[0]);
    calculate_and_print(*veem, prepared[0]);

    const DemoAtomParameterMethod demo_method;
    const auto demo_options = methods::make_default_options(demo_method.option_schema());

    print_result("\ndemo method prerequisites for collection",
                 methods::check_method_prerequisites(demo_method, prepared, demo_options));

    const auto water_parameters = make_water_parameters();

    print_parameter_result(
        "demo parameter prerequisites for first molecule",
        methods::check_parameter_prerequisites(demo_method, {.prepared_molecule = prepared[0],
                                                             .parameter_set = water_parameters,
                                                             .classification_options = {}}));

    print_collection_parameter_result(
        "demo parameter prerequisites for whole collection",
        methods::check_parameter_prerequisites(demo_method, prepared, water_parameters, {}));

    return 0;
}