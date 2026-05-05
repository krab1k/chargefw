#include "methods/builtin/delre.h"

#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/LU>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

[[nodiscard]] auto matches_atomic_number(const core::Atom& atom,
                                         const parameters::AtomParameterKey& key) -> bool {
    return key.atomic_number == 0 || key.atomic_number == atom.atomic_number();
}

[[nodiscard]] auto bond_parameters_are_reversed(const core::Molecule& molecule,
                                                const parameters::ParameterView& parameters,
                                                const std::size_t bond_index) -> bool {
    const auto& bond = molecule.bond(bond_index);
    const auto& first_atom = molecule.atom(bond.first_atom_index());
    const auto& second_atom = molecule.atom(bond.second_atom_index());

    const auto parameter_entry_index =
        parameters.classification().bond().parameter_entry_index(bond_index);
    const auto& key = parameters.parameter_set().bond().entry(parameter_entry_index).key;

    const auto forward = matches_atomic_number(first_atom, key.first_atom) &&
                         matches_atomic_number(second_atom, key.second_atom);

    if (forward) {
        return false;
    }

    const auto reverse = matches_atomic_number(first_atom, key.second_atom) &&
                         matches_atomic_number(second_atom, key.first_atom);

    if (reverse) {
        return true;
    }

    throw std::logic_error{"DelRe bond parameter classification is inconsistent"};
}

} // namespace

auto DelReMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto delta = parameters.atom("delta");
    const auto eps = parameters.bond("eps");
    const auto gamma_a = parameters.bond("gammaA");
    const auto gamma_b = parameters.bond("gammaB");

    const auto n = static_cast<Eigen::Index>(atom_count);

    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto eigen_index = static_cast<Eigen::Index>(atom_index);

        rhs(eigen_index) = -delta[atom_index];
        matrix(eigen_index, eigen_index) = -1.0;
    }

    for (std::size_t bond_index = 0; bond_index < molecule.bond_count(); ++bond_index) {
        const auto& bond = molecule.bond(bond_index);
        const auto first_index = static_cast<Eigen::Index>(bond.first_atom_index());
        const auto second_index = static_cast<Eigen::Index>(bond.second_atom_index());

        const auto reversed = bond_parameters_are_reversed(molecule, parameters, bond_index);

        matrix(first_index, second_index) = reversed ? gamma_b[bond_index] : gamma_a[bond_index];
        matrix(second_index, first_index) = reversed ? gamma_a[bond_index] : gamma_b[bond_index];
    }

    const Eigen::VectorXd solution = matrix.partialPivLu().solve(rhs);

    auto values = std::vector(atom_count, 0.0);

    for (std::size_t bond_index = 0; bond_index < molecule.bond_count(); ++bond_index) {
        const auto& bond = molecule.bond(bond_index);
        const auto first_atom_index = bond.first_atom_index();
        const auto second_atom_index = bond.second_atom_index();

        const auto denominator = 2.0 * eps[bond_index];

        if (denominator == 0.0) {
            throw std::logic_error{"DelRe bond parameter 'eps' must be non-zero"};
        }

        const auto first_index = static_cast<Eigen::Index>(first_atom_index);
        const auto second_index = static_cast<Eigen::Index>(second_atom_index);

        const auto bond_charge = (solution(first_index) - solution(second_index)) / denominator;

        values[first_atom_index] -= bond_charge;
        values[second_atom_index] += bond_charge;
    }

    return charges::AtomicCharges{std::move(values)};
}

} // namespace chargefw::methods::builtin