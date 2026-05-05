#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <exception>
#include <iostream>
#include <string>
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

    std::vector positions_1{core::Position{.x = 0.0000, .y = 0.0000, .z = 0.0000},
                            core::Position{.x = 0.9572, .y = 0.0000, .z = 0.0000},
                            core::Position{.x = -0.2390, .y = 0.9270, .z = 0.0000}};

    std::vector positions_2{core::Position{.x = 0.0000, .y = 0.0000, .z = 0.0000},
                            core::Position{.x = 1.1000, .y = 0.0000, .z = 0.0000},
                            core::Position{.x = -0.3000, .y = 1.0500, .z = 0.0000}};

    std::vector conformers{core::Conformer{std::move(positions_1), "model-1"},
                           core::Conformer{std::move(positions_2), "model-2"}};

    return core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers), "water"};
}

auto make_collection() -> core::MoleculeCollection {
    std::vector molecules{make_water()};

    return core::MoleculeCollection{std::move(molecules), "demo"};
}

auto method_pointers(const methods::MethodRegistry& registry)
    -> std::vector<const methods::Method*> {
    std::vector<const methods::Method*> result;
    result.reserve(registry.methods().size());

    for (const auto& method : registry.methods()) {
        result.push_back(method.get());
    }

    return result;
}

auto print_charge_set(const charges::ChargeSet& charge_set) -> void {
    std::cout << "method: " << charge_set.method_id();

    if (charge_set.parameter_set_id().has_value()) {
        std::cout << "  parameters: " << *charge_set.parameter_set_id();
    }

    std::cout << '\n';

    for (const auto& assignment : charge_set.assignments()) {
        std::cout << "  molecule " << assignment.target.molecule_index;

        if (assignment.target.conformer_index.has_value()) {
            std::cout << " conformer " << *assignment.target.conformer_index;
        }

        std::cout << " charges:";

        for (const auto charge : assignment.charges.values()) {
            std::cout << ' ' << charge;
        }

        std::cout << "  total=" << assignment.charges.total() << '\n';
    }
}

auto print_rejections(const methods::ApplicabilityResult& applicability,
                      const std::vector<const methods::Method*>& candidate_methods,
                      const std::vector<parameters::ParameterSet>& parameter_sets) -> void {
    if (applicability.rejected.empty()) {
        return;
    }

    std::cout << "\nRejected candidates:\n";

    for (const auto& rejected : applicability.rejected) {
        const auto* method = candidate_methods[rejected.method_index];

        std::cout << "  method: " << method->id();

        if (rejected.parameter_set_index.has_value()) {
            const auto& parameter_set = parameter_sets[*rejected.parameter_set_index];
            std::cout << "  parameters: " << parameter_set.id();
        }

        std::cout << '\n';

        for (const auto& issue : rejected.issues) {
            std::cout << "    - " << issue.message << '\n';
        }
    }
}

} // namespace

auto main() -> int {
    try {
        const auto collection = make_collection();
        const features::PreparedMoleculeCollection prepared_collection{collection};

        const auto parameter_sets = parameters::load_default_parameter_sets();

        const auto& registry = methods::method_registry();
        const auto candidates = method_pointers(registry);

        const auto applicability =
            methods::find_applicable_methods(prepared_collection, candidates, parameter_sets);

        std::cout << "Loaded methods: " << candidates.size() << '\n';
        std::cout << "Loaded parameter sets: " << parameter_sets.size() << '\n';
        std::cout << "Applicable candidates: " << applicability.applicable.size() << "\n\n";

        for (const auto& candidate : applicability.applicable) {
            try {
                const auto charge_set = methods::calculate_charges(candidate, prepared_collection);
                print_charge_set(charge_set);
                std::cout << '\n';
            } catch (const std::exception& error) {
                std::cerr << "Calculation failed for method '" << candidate.method->id() << "'";

                if (candidate.parameter_set != nullptr) {
                    std::cerr << " with parameter set '" << candidate.parameter_set->id() << "'";
                }

                std::cerr << ": " << error.what() << '\n';
            }
        }

        print_rejections(applicability, candidates, parameter_sets);

        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Fatal error: " << error.what() << '\n';
        return 1;
    }
}