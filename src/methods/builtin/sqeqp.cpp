#include "methods/builtin/sqeqp.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/LU>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

[[nodiscard]] auto interaction(const double distance, const double width_i, const double width_j)
    -> double {
    const auto width_sum = 2.0 * width_i * width_i + 2.0 * width_j * width_j;

    if (width_sum == 0.0) {
        return 1.0 / distance;
    }

    return std::erf(distance / std::sqrt(width_sum)) / distance;
}

} // namespace

auto SQEqpMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& geometry = input.geometry();
    const auto& parameters = input.parameters();

    const auto atom_count = molecule.atom_count();
    const auto bond_count = molecule.bond_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto electronegativity = parameters.atom("electronegativity");
    const auto hardness = parameters.atom("hardness");
    const auto width = parameters.atom("width");
    const auto initial_charge_parameter = parameters.atom("q0");

    const auto n = static_cast<Eigen::Index>(atom_count);
    Eigen::VectorXd initial_charges = Eigen::VectorXd::Zero(n);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        initial_charges(static_cast<Eigen::Index>(atom_index)) =
            initial_charge_parameter[atom_index];
    }

    initial_charges.array() -= (initial_charges.sum() - core::total_formal_charge(molecule)) /
                               static_cast<double>(atom_count);

    if (bond_count == 0) {
        return charges::AtomicCharges{
            std::vector<double>{initial_charges.data(), initial_charges.data() + atom_count}};
    }

    const auto kappa = parameters.bond("kappa");
    const auto m = static_cast<Eigen::Index>(bond_count);

    Eigen::MatrixXd charge_matrix = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd charge_rhs = Eigen::VectorXd::Zero(n);
    Eigen::MatrixXd transfer_matrix = Eigen::MatrixXd::Zero(m, n);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);
        charge_matrix(i, i) = hardness[atom_index];
        charge_rhs(i) = -electronegativity[atom_index];

        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto value = interaction(geometry.distance(atom_index, other_atom_index),
                                           width[atom_index], width[other_atom_index]);
            charge_matrix(i, j) = value;
            charge_matrix(j, i) = value;
        }
    }

    for (std::size_t bond_index = 0; bond_index < bond_count; ++bond_index) {
        const auto& bond = molecule.bond(bond_index);
        const auto row = static_cast<Eigen::Index>(bond_index);
        transfer_matrix(row, static_cast<Eigen::Index>(bond.first_atom_index())) = 1.0;
        transfer_matrix(row, static_cast<Eigen::Index>(bond.second_atom_index())) = -1.0;
    }

    charge_rhs -= charge_matrix * initial_charges;
    charge_rhs += charge_matrix.diagonal().cwiseProduct(initial_charges);

    Eigen::MatrixXd split_matrix = transfer_matrix * charge_matrix * transfer_matrix.transpose();

    for (std::size_t bond_index = 0; bond_index < bond_count; ++bond_index) {
        split_matrix(static_cast<Eigen::Index>(bond_index),
                     static_cast<Eigen::Index>(bond_index)) += kappa[bond_index];
    }

    const Eigen::VectorXd split_charge =
        split_matrix.partialPivLu().solve(transfer_matrix * charge_rhs);
    const Eigen::VectorXd charges = transfer_matrix.transpose() * split_charge + initial_charges;

    return charges::AtomicCharges{
        std::vector<double>{charges.data(), charges.data() + charges.size()}};
}

} // namespace chargefw::methods::builtin
