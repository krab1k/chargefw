#include "methods/builtin/kcm.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/LU>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {

auto KCMMethod::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                        PrerequisiteResult& result) const -> void {
    if (core::total_formal_charge(input.prepared_molecule.molecule()) != 0) {
        result.add(PrerequisiteIssue{
            .kind = PrerequisiteIssueKind::unsupported_molecule,
            .message =
                "KCM supports only neutral molecules because counter-ion graph augmentation is not "
                "implemented"});
    }
}

auto KCMMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();

    const auto n = static_cast<Eigen::Index>(molecule.atom_count());

    if (n == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto electronegativity = parameters.atom("electronegativity");
    const auto hardness = parameters.atom("hardness");

    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n, n);
    Eigen::VectorXd chi0 = Eigen::VectorXd::Zero(n);

    for (Eigen::Index i = 0; i < n; ++i) {
        chi0(i) = electronegativity[static_cast<std::size_t>(i)];
    }

    for (const auto& bond : molecule.bonds()) {
        const auto i = bond.first_atom_index();
        const auto j = bond.second_atom_index();

        const auto eta = hardness[i] + hardness[j];

        if (eta == 0.0) {
            throw std::logic_error{"KCM bond hardness denominator must be non-zero"};
        }

        const auto w = 1.0 / eta;

        const auto ei = static_cast<Eigen::Index>(i);
        const auto ej = static_cast<Eigen::Index>(j);

        A(ei, ei) += w;
        A(ej, ej) += w;
        A(ei, ej) -= w;
        A(ej, ei) -= w;
    }

    Eigen::VectorXd q = A.partialPivLu().solve(chi0);
    q -= chi0;

    return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
}

} // namespace chargefw::methods::builtin
