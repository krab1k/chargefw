#include "methods/builtin/sfkeem.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/LU>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

[[nodiscard]] auto interaction(const double sigma, const double distance, const double b_i,
                               const double b_j) -> double {

    if (b_i <= 0.0 || b_j <= 0.0) {
        throw std::logic_error{"SFKEEM atom parameter 'B' must be positive"};
    }

    return 2.0 * std::sqrt(b_i * b_j) / std::cosh(sigma * distance);
}

} // namespace

auto SFKEEMMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();
    const auto& geometry = input.geometry();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto sigma = parameters.common("sigma");

    if (sigma <= 0.0) {
        throw std::logic_error{"SFKEEM common parameter 'sigma' must be positive"};
    }

    const auto a = parameters.atom("A");
    const auto b_parameter = parameters.atom("B");

    const auto n = static_cast<Eigen::Index>(atom_count);

    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + 1);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);

        matrix(i, i) = 2.0 * b_parameter[atom_index];
        rhs(i) = -a[atom_index];

        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto distance = geometry.distance(atom_index, other_atom_index);

            const auto value = interaction(sigma, distance, b_parameter[atom_index],
                                           b_parameter[other_atom_index]);

            matrix(i, j) = value;
            matrix(j, i) = value;
        }
    }

    matrix.row(n) = Eigen::VectorXd::Constant(n + 1, 1.0);
    matrix.col(n) = Eigen::VectorXd::Constant(n + 1, 1.0);
    matrix(n, n) = 0.0;
    rhs(n) = input.target_charge();

    const Eigen::VectorXd q = matrix.partialPivLu().solve(rhs).head(n);

    return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
}

} // namespace chargefw::methods::builtin
