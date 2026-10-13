#include "methods/builtin/denr.h"

#include "methods/builtin/prerequisite_helpers.h"

#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/LU>

#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {

auto DENRMethod::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                         PrerequisiteResult& result) const -> void {
    detail::add_component_neutrality_prerequisite_issue(input, result, "DENR");
}

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
    const auto iterations = options.get<int>("iters");

    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd A = I + dt * L * eta;

    // Factorize in place to avoid holding a second dense copy of A.
    const Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXd>> lu{A};
    const Eigen::VectorXd tmp = dt * L * chi;

    for (auto i = 0; i < iterations; ++i) {
        q = lu.solve(q - tmp);
    }

    return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
}

auto DENRMethod::estimate(const TargetSize& target, const MethodOptions& options) const
    -> ResourceEstimate {
    const auto n = target.atoms;
    const auto iterations = options.get<int>("iters");
    // eta, L, I, and A are dense n x n; L * eta, LU, and one solve per iteration.
    return {.memory_bytes = 32.0 * n * n,
            .operations = 2.0 * n * n * n + 2.0 / 3.0 * n * n * n + iterations * 2.0 * n * n};
}

} // namespace chargefw::methods::builtin
