#include "methods/builtin/kcm.h"

#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/LU>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {

auto KCMMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();

    const auto n = static_cast<Eigen::Index>(molecule.atom_count());
    const auto m = static_cast<Eigen::Index>(molecule.bond_count());

    if (n == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto electronegativity = parameters.atom("electronegativity");
    const auto hardness = parameters.atom("hardness");

    Eigen::MatrixXd W = Eigen::MatrixXd::Zero(m, m);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m, n);
    Eigen::VectorXd chi0 = Eigen::VectorXd::Zero(n);

    for (Eigen::Index i = 0; i < n; ++i) {
        chi0(i) = electronegativity[static_cast<std::size_t>(i)];
    }

    for (Eigen::Index i = 0; i < m; ++i) {
        const auto bond_index = static_cast<std::size_t>(i);
        const auto& bond = molecule.bond(bond_index);

        const auto first = bond.first_atom_index();
        const auto second = bond.second_atom_index();

        const auto eta = hardness[first] + hardness[second];

        if (eta == 0.0) {
            throw std::logic_error{"KCM bond hardness denominator must be non-zero"};
        }

        W(i, i) = 1.0 / eta;

        B(i, static_cast<Eigen::Index>(first)) = 1.0;
        B(i, static_cast<Eigen::Index>(second)) = -1.0;
    }

    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
    const Eigen::MatrixXd A = B.transpose() * W * B + I;

    Eigen::VectorXd q = A.partialPivLu().solve(chi0);
    q -= chi0;

    return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
}

} // namespace chargefw::methods::builtin