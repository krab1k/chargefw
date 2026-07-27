#include "methods/builtin/mpeoe.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

[[nodiscard]] auto is_hydrogen(const core::Atom& atom) noexcept -> bool {
    return atom.atomic_number() == 1;
}

[[nodiscard]] auto atom_denominator(const parameters::AtomParameterAccessor& parameter_a,
                                    const parameters::AtomParameterAccessor& parameter_b,
                                    const std::size_t atom_index) -> double {
    return parameter_a[atom_index] + parameter_b[atom_index];
}

} // namespace

auto MPEOEMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto iterations = input.method_options().get<int>("iters");

    if (iterations <= 0) {
        throw std::invalid_argument{"MPEOE option 'iters' must be positive"};
    }

    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();

    const auto parameter_a = parameters.atom("A");
    const auto parameter_b = parameters.atom("B");
    const auto parameter_f = parameters.bond("f");
    const auto hydrogen_plus = parameters.common("Hplus");

    const auto atom_count = molecule.atom_count();

    std::vector charges(atom_count, 0.0);
    std::vector electronegativities(atom_count, 0.0);

    for (int iteration = 0; iteration < iterations; ++iteration) {
        const auto alpha = iteration + 1;
        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            electronegativities[atom_index] =
                parameter_b[atom_index] * charges[atom_index] + parameter_a[atom_index];
        }

        for (std::size_t bond_index = 0; bond_index < molecule.bond_count(); ++bond_index) {
            const auto& bond = molecule.bond(bond_index);

            auto lower_atom_index = bond.first_atom_index();
            auto higher_atom_index = bond.second_atom_index();

            auto lower_electronegativity = electronegativities[lower_atom_index];
            auto higher_electronegativity = electronegativities[higher_atom_index];

            if (lower_electronegativity > higher_electronegativity) {
                std::swap(lower_atom_index, higher_atom_index);
                std::swap(lower_electronegativity, higher_electronegativity);
            }

            const auto& lower_atom = molecule.atom(lower_atom_index);

            const auto denominator =
                is_hydrogen(lower_atom)
                    ? hydrogen_plus
                    : atom_denominator(parameter_a, parameter_b, lower_atom_index);

            const auto difference = std::pow(parameter_f[bond_index], alpha) *
                                    (higher_electronegativity - lower_electronegativity) /
                                    denominator;

            charges[lower_atom_index] += difference;
            charges[higher_atom_index] -= difference;
        }
    }

    return charges::AtomicCharges{std::move(charges)};
}

} // namespace chargefw::methods::builtin