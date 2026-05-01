#include "methods/builtin/peoe.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/parameters/parameter_view.h>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

[[nodiscard]] auto is_hydrogen(const core::Atom& atom) noexcept -> bool {
    return atom.atomic_number() == 1;
}

[[nodiscard]] auto atom_denominator(const parameters::AtomParameterAccessor& parameter_a,
                                    const parameters::AtomParameterAccessor& parameter_b,
                                    const parameters::AtomParameterAccessor& parameter_c,
                                    const std::size_t atom_index) -> double {
    return parameter_a[atom_index] + parameter_b[atom_index] + parameter_c[atom_index];
}

} // namespace

auto PEOEMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto iterations = input.method_options().get<int>("iters");

    if (iterations <= 0) {
        throw std::invalid_argument{"PEOE option 'iters' must be positive"};
    }

    const auto& molecule = input.molecule();
    const auto& parameter_view = input.parameters();

    const auto parameter_a = parameter_view.atom("A");
    const auto parameter_b = parameter_view.atom("B");
    const auto parameter_c = parameter_view.atom("C");
    const auto damping_hydrogen = parameter_view.common("dampH");

    const auto atom_count = molecule.atom_count();

    std::vector<double> charges(atom_count, 0.0);
    std::vector<double> electronegativities(atom_count, 0.0);

    for (int alpha = 1; alpha < iterations; ++alpha) {
        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            const auto charge = charges[atom_index];

            electronegativities[atom_index] = parameter_c[atom_index] * charge * charge +
                                              parameter_b[atom_index] * charge +
                                              parameter_a[atom_index];
        }

        for (const auto& bond : molecule.bonds()) {
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
                    ? damping_hydrogen
                    : atom_denominator(parameter_a, parameter_b, parameter_c, lower_atom_index);

            const auto difference = std::pow(0.5, alpha) *
                                    (higher_electronegativity - lower_electronegativity) /
                                    denominator;

            charges[lower_atom_index] += difference;
            charges[higher_atom_index] -= difference;
        }
    }

    return charges::AtomicCharges{std::move(charges)};
}

} // namespace chargefw::methods::builtin