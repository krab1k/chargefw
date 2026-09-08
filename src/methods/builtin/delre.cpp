#include "methods/builtin/delre.h"

#include "features/topology_helpers.h"

#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/LU>

#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {

auto DelReMethod::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                          PrerequisiteResult& result) const
    -> void {
    const auto& molecule = input.prepared_molecule.molecule();
    const auto components =
        features::connected_components(input.prepared_molecule.topology().adjacency());

    for (const auto& component : components) {
        const auto formal_charge =
            std::accumulate(component.begin(), component.end(), 0,
                            [&molecule](const int sum, const std::size_t atom_index) {
                                return sum + molecule.atom(atom_index).formal_charge();
                            });
        if (formal_charge != 0) {
            result.add(PrerequisiteIssue{
                .kind = PrerequisiteIssueKind::unsupported_molecule,
                .message =
                    "DelRe supports only molecules with neutral connected components because its "
                    "bond-charge construction conserves zero total charge within each component"});
            return;
        }
    }
}

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

        const auto parameter_entry_index =
            parameters.classification().bond().parameter_entry_index(bond_index);
        const auto& key = parameters.parameter_set().bond().entry(parameter_entry_index).key;
        // DelRe_original has plain element keys and symmetric homonuclear gamma values.
        const auto first_is_parameter_a =
            molecule.atom(bond.first_atom_index()).atomic_number() == key.first_atom.atomic_number;

        matrix(first_index, second_index) =
            first_is_parameter_a ? gamma_a[bond_index] : gamma_b[bond_index];
        matrix(second_index, first_index) =
            first_is_parameter_a ? gamma_b[bond_index] : gamma_a[bond_index];
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
