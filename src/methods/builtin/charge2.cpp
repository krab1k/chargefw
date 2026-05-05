#include "methods/builtin/charge2.h"

#include "features/topology_helpers.h"

#include <chargefw/core/periodic_table.h>
#include <chargefw/parameters/models/parameter_view.h>
#include <math.h>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

constexpr auto hydrogen_atomic_number = 1;
constexpr auto max_bond_distance = std::size_t{3};

struct Charge2Shells {
    std::vector<std::size_t> one_bond;
    std::vector<std::size_t> two_bonds;
    std::vector<std::size_t> three_bonds;
};

[[nodiscard]] auto precompute_charge2_shells(const std::vector<std::vector<std::size_t>>& adjacency)
    -> std::vector<Charge2Shells> {
    auto shells_by_atom = std::vector<Charge2Shells>(adjacency.size());

    for (std::size_t atom_index = 0; atom_index < adjacency.size(); ++atom_index) {
        const auto distances =
            features::breadth_first_distances(adjacency, atom_index, max_bond_distance);

        auto& shells = shells_by_atom[atom_index];

        for (std::size_t bonded_index = 0; bonded_index < distances.size(); ++bonded_index) {
            if (distances[bonded_index] == 1) {
                shells.one_bond.push_back(bonded_index);
            } else if (distances[bonded_index] == 2) {
                shells.two_bonds.push_back(bonded_index);
            } else if (distances[bonded_index] == 3) {
                shells.three_bonds.push_back(bonded_index);
            }
        }
    }

    return shells_by_atom;
}

[[nodiscard]] auto hydrogen_electronegativity(const parameters::ParameterView& parameters)
    -> double {
    const auto& atom_parameters = parameters.parameter_set().atom();

    for (std::size_t parameter_index = 0; parameter_index < atom_parameters.size();
         ++parameter_index) {
        const auto& entry = atom_parameters[parameter_index];

        if (entry.key.atomic_number == hydrogen_atomic_number &&
            atom_parameters.contains(parameter_index, "chi")) {
            return atom_parameters.parameter(parameter_index, "chi");
        }
    }

    throw std::logic_error{"Charge2 requires hydrogen atom parameter 'chi'"};
}

} // namespace

auto Charge2Method::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto iterations = input.method_options().get<int>("iters");

    if (iterations <= 0) {
        throw std::invalid_argument{"Charge2 option 'iters' must be positive"};
    }

    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();
    const auto& table = core::periodic_table();

    const auto chi = parameters.atom("chi");
    const auto p0 = parameters.atom("P0");
    const auto q0 = parameters.atom("q0");

    const auto a1 = parameters.common("a1");
    const auto a2 = parameters.common("a2");
    const auto a3 = parameters.common("a3");
    const auto b = parameters.common("b");
    const auto c = parameters.common("c");
    const auto alpha = parameters.common("alpha");

    const auto chi_hydrogen = hydrogen_electronegativity(parameters);

    const auto atom_count = molecule.atom_count();
    const auto shells_by_atom = precompute_charge2_shells(input.topology().adjacency());

    auto q = std::vector(atom_count, 0.0);

    for (auto iteration = 0; iteration < iterations; ++iteration) {
        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            const auto& atom = molecule.atom(atom_index);
            const auto& [one_bond_apart, two_bonds_apart, three_bonds_apart] =
                shells_by_atom[atom_index];

            auto alpha_charge = 0.0;

            for (const auto bonded_index : one_bond_apart) {
                const auto& bonded = molecule.atom(bonded_index);

                double a = 0.0;
                if (table.element(atom.atomic_number()).period == 2 &&
                    table.element(bonded.atomic_number()).period == 2) {
                    a = a1;
                } else if (atom.atomic_number() == hydrogen_atomic_number ||
                           bonded.atomic_number() == hydrogen_atomic_number) {
                    a = a2;
                } else {
                    a = a3;
                }

                const auto bonded_electronegativity = chi[bonded_index];
                const auto atom_electronegativity = chi[atom_index];

                alpha_charge += (bonded_electronegativity - atom_electronegativity) / a;
            }

            const auto polarizability =
                p0[atom_index] * (1.0 + alpha * (q0[atom_index] - q[atom_index]));

            auto beta_charge = 0.0;

            for (const auto bonded_index : two_bonds_apart) {
                beta_charge += (chi[bonded_index] - chi_hydrogen) * polarizability / b;
            }

            auto gamma_charge = 0.0;

            for (const auto bonded_index : three_bonds_apart) {
                gamma_charge += (chi[bonded_index] - chi_hydrogen) * polarizability / b / c;
            }

            q[atom_index] = alpha_charge + beta_charge + gamma_charge;
        }
    }

    return charges::AtomicCharges{std::move(q)};
}

} // namespace chargefw::methods::builtin