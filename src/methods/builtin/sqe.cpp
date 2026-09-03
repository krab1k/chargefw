#include "methods/builtin/sqe.h"

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

auto sqe_core::calculate(const CalculationInput& input,
                         const std::span<const double> initial_charge_values)
    -> std::vector<double> {
    const auto& molecule = input.molecule();
    const auto& geometry = input.geometry();
    const auto& parameters = input.parameters();

    const auto atom_count = molecule.atom_count();
    const auto bond_count = molecule.bond_count();

    if (atom_count == 0) {
        return {};
    }

    if (!initial_charge_values.empty() && initial_charge_values.size() != atom_count) {
        throw std::logic_error{"SQE-family initial charge count must match atom count"};
    }

    if (bond_count == 0) {
        if (initial_charge_values.empty()) {
            return std::vector<double>(atom_count, 0.0);
        }

        return {initial_charge_values.begin(), initial_charge_values.end()};
    }

    const auto electronegativity = parameters.atom("electronegativity");
    const auto hardness = parameters.atom("hardness");
    const auto width = parameters.atom("width");
    const auto kappa = parameters.bond("kappa");

    const auto n = static_cast<Eigen::Index>(atom_count);
    const auto m = static_cast<Eigen::Index>(bond_count);

    Eigen::MatrixXd charge_matrix = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd charge_rhs = Eigen::VectorXd::Zero(n);
    Eigen::MatrixXd transfer_matrix = Eigen::MatrixXd::Zero(m, n);
    Eigen::VectorXd initial_charges = Eigen::VectorXd::Zero(n);

    if (!initial_charge_values.empty()) {
        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            initial_charges(static_cast<Eigen::Index>(atom_index)) =
                initial_charge_values[atom_index];
        }
    }

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

    if (!initial_charge_values.empty()) {
        charge_rhs -= charge_matrix * initial_charges;
        charge_rhs += charge_matrix.diagonal().cwiseProduct(initial_charges);
    }

    Eigen::MatrixXd split_matrix = transfer_matrix * charge_matrix * transfer_matrix.transpose();

    for (std::size_t bond_index = 0; bond_index < bond_count; ++bond_index) {
        split_matrix(static_cast<Eigen::Index>(bond_index),
                     static_cast<Eigen::Index>(bond_index)) += kappa[bond_index];
    }

    const Eigen::VectorXd split_charge =
        split_matrix.partialPivLu().solve(transfer_matrix * charge_rhs);
    Eigen::VectorXd charges = transfer_matrix.transpose() * split_charge;

    if (!initial_charge_values.empty()) {
        charges += initial_charges;
    }

    return {charges.data(), charges.data() + charges.size()};
}

auto SQEMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    return charges::AtomicCharges{sqe_core::calculate(input)};
}

} // namespace chargefw::methods::builtin
