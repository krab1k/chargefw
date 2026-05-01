#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/parameter_set_io.h>

#include <iostream>
#include <utility>
#include <vector>

namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

auto make_hf() -> core::Molecule {
    std::vector atoms{core::Atom{1, 0, "H"}, core::Atom{9, 0, "F"}};

    std::vector bonds{core::Bond{0, 1, core::BondOrder::SINGLE}};

    return core::Molecule{std::move(atoms), std::move(bonds), {}, "hf"};
}

auto make_collection() -> core::MoleculeCollection {
    std::vector molecules{make_hf()};

    return core::MoleculeCollection{std::move(molecules), "demo"};
}

auto print_charge_set(const charges::ChargeSet& charge_set) -> void {
    std::cout << "method: " << charge_set.method_id() << '\n';

    if (charge_set.parameter_set_id().has_value()) {
        std::cout << "parameters: " << *charge_set.parameter_set_id() << '\n';
    }

    for (const auto& assignment : charge_set.assignments()) {
        std::cout << "molecule " << assignment.target.molecule_index << " charges:";

        for (const auto charge : assignment.charges.values()) {
            std::cout << ' ' << charge;
        }

        std::cout << "  total=" << assignment.charges.total() << '\n';
    }
}

} // namespace

auto main() -> int {
    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared{collection};

    const auto parameter_sets = parameters::load_default_parameter_sets();

    if (parameter_sets.empty()) {
        std::cerr << "No default parameter sets found.\n";
        return 1;
    }

    const auto& registry = methods::method_registry();
    const auto* peoe = registry.find("peoe");

    if (peoe == nullptr) {
        std::cerr << "PEOE method is not registered.\n";
        return 1;
    }

    const std::vector candidate_methods{peoe};

    const auto applicability =
        methods::find_applicable_methods(prepared, candidate_methods, parameter_sets);

    if (applicability.applicable.empty()) {
        std::cerr << "No applicable method/parameter-set pair.\n";

        for (const auto& rejected : applicability.rejected) {
            for (const auto& issue : rejected.issues) {
                std::cerr << "  - " << issue.message << '\n';
            }
        }

        return 1;
    }

    const auto charge_set = methods::calculate_charges(applicability.applicable.front(), prepared);

    print_charge_set(charge_set);

    return 0;
}