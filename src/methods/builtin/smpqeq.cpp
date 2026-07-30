#include "methods/builtin/smpqeq.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/LU>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

constexpr auto iteration_count = 5;

[[nodiscard]] auto interaction(const double distance, const double second_i, const double second_j)
    -> double {

    if (second_i <= 0.0 || second_j <= 0.0) {
        throw std::logic_error{"SMP/QEq atom parameter 'second' must be positive"};
    }

    const auto gamma = 2.0 * std::sqrt(second_i * second_j);

    return 1.0 / std::cbrt(1.0 / std::pow(gamma, 3.0) + std::pow(distance, 3.0));
}

[[nodiscard]] auto diagonal_term(const double previous_charge, const double second,
                                 const double third, const double fourth) -> double {
    return 2.0 * (second + third * previous_charge + fourth * previous_charge * previous_charge);
}

} // namespace

auto SMPQEqMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();
    const auto& geometry = input.geometry();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto first = parameters.atom("first");
    const auto second = parameters.atom("second");
    const auto third = parameters.atom("third");
    const auto fourth = parameters.atom("fourth");

    const auto n = static_cast<Eigen::Index>(atom_count);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);

    for (auto iteration = 0; iteration < iteration_count; ++iteration) {
        A.setZero();

        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            const auto i = static_cast<Eigen::Index>(atom_index);
            const auto previous_charge = b(i);

            A(i, i) = diagonal_term(previous_charge, second[atom_index], third[atom_index],
                                    fourth[atom_index]);
            b(i) = -first[atom_index];

            for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
                 ++other_atom_index) {
                const auto j = static_cast<Eigen::Index>(other_atom_index);
                const auto distance = geometry.distance(atom_index, other_atom_index);

                const auto value =
                    interaction(distance, second[atom_index], second[other_atom_index]);

                A(i, j) = value;
                A(j, i) = value;
            }
        }

        A.row(n) = Eigen::VectorXd::Constant(n + 1, 1.0);
        A.col(n) = Eigen::VectorXd::Constant(n + 1, 1.0);
        A(n, n) = 0.0;
        b(n) = input.target_charge();

        b = A.partialPivLu().solve(b);
    }

    return charges::AtomicCharges{std::vector<double>{b.data(), b.data() + n}};
}

} // namespace chargefw::methods::builtin
