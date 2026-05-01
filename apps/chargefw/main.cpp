#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/atom_parameters.h>
#include <chargefw/parameters/parameter_key.h>
#include <chargefw/parameters/parameter_set.h>
#include <chargefw/parameters/parameter_set_metadata.h>

#include <cstddef>
#include <iostream>
#include <optional>
#include <span>
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

auto make_demo_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "demo-parameters",
                                         .method_id = "demo-atom-parameter-method",
                                         .name = "Demo parameters"},
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
        const auto value = input.parameters().atom("value");

        std::vector<double> values;
        values.reserve(molecule.atom_count());

        for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
            values.push_back(value[atom_index]);
        }

        return charges::AtomicCharges{std::move(values)};
    }
};

auto print_applicable(const methods::ApplicabilityResult& result) -> void {
    std::cout << "Applicable candidates: " << result.applicable.size() << '\n';

    for (const auto& candidate : result.applicable) {
        std::cout << "  - method: " << candidate.method->id();

        if (candidate.parameter_set != nullptr) {
            std::cout << ", parameters: " << candidate.parameter_set->id();
            std::cout << ", classifications: " << candidate.classifications.size();
        } else {
            std::cout << ", no parameters";
        }

        std::cout << '\n';
    }
}

auto select_demo_candidate(const methods::ApplicabilityResult& result)
    -> const methods::ApplicableMethod* {
    for (const auto& candidate : result.applicable) {
        if (candidate.uses_parameters()) {
            return &candidate;
        }
    }

    if (!result.applicable.empty()) {
        return &result.applicable.front();
    }

    return nullptr;
}

auto print_charge_set(const charges::ChargeSet& charge_set) -> void {
    std::cout << "\nCharge set:\n";
    std::cout << "  method: " << charge_set.method_id() << '\n';

    if (charge_set.parameter_set_id().has_value()) {
        std::cout << "  parameters: " << *charge_set.parameter_set_id() << '\n';
    }

    for (const auto& assignment : charge_set.assignments()) {
        std::cout << "  molecule " << assignment.target.molecule_index << " charges:";

        for (const auto charge : assignment.charges.values()) {
            std::cout << ' ' << charge;
        }

        std::cout << "  total=" << assignment.charges.total() << '\n';
    }
}

} // namespace

auto main() -> int {
    const auto collection = make_demo_collection();
    const features::PreparedMoleculeCollection prepared{collection};

    const auto& registry = methods::method_registry();

    std::vector<const methods::Method*> candidate_methods;

    for (const auto& method : registry.methods()) {
        candidate_methods.push_back(method.get());
    }

    const DemoAtomParameterMethod demo_method;
    candidate_methods.push_back(&demo_method);

    const std::vector parameter_sets{make_demo_parameters()};

    const auto applicability =
        methods::find_applicable_methods(prepared, candidate_methods, parameter_sets);

    print_applicable(applicability);

    const auto* selected = select_demo_candidate(applicability);

    if (selected == nullptr) {
        std::cerr << "\nNo method or method/parameter-set pair supports the whole collection.\n";
        return 1;
    }

    const auto charge_set = methods::calculate_charges(*selected, prepared);

    print_charge_set(charge_set);

    return 0;
}