#include "methods/builtin/denr.h"

#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/LU>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {

auto DENRMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();
    const auto& options = input.method_options();

    const auto n = static_cast<Eigen::Index>(molecule.atom_count());

    if (n == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto electronegativity = parameters.atom("electronegativity");
    const auto hardness = parameters.atom("hardness");

    Eigen::MatrixXd eta = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd chi = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd q = Eigen::VectorXd::Zero(n);

    for (Eigen::Index i = 0; i < n; ++i) {
        const auto atom_index = static_cast<std::size_t>(i);

        chi(i) = electronegativity[atom_index];
        eta(i, i) = hardness[atom_index];
    }

    for (const auto& bond : molecule.bonds()) {
        const auto i = static_cast<Eigen::Index>(bond.first_atom_index());
        const auto j = static_cast<Eigen::Index>(bond.second_atom_index());

        L(i, i) += 1.0;
        L(j, j) += 1.0;
        L(i, j) -= 1.0;
        L(j, i) -= 1.0;
    }

    const auto dt = options.get<double>("step");
    const auto iterations = options.get<int>("iterations");

    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
    const Eigen::MatrixXd A = I + dt * L * eta;

    const auto lu = A.partialPivLu();
    const Eigen::VectorXd tmp = dt * L * chi;

    for (auto i = 0; i < iterations; ++i) {
        q = lu.solve(q - tmp);
    }

    return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
}

} // namespace chargefw::methods::builtin
