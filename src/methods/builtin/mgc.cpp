#include "methods/builtin/mgc.h"

#include <chargefw/core/periodic_table.h>

#include <Eigen/LU>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

[[nodiscard]] auto order(const core::BondOrder bond_order) -> double {
    switch (bond_order) {
    case core::BondOrder::UNKNOWN:
    case core::BondOrder::SINGLE:
        return 1.0;
    case core::BondOrder::DOUBLE:
        return 2.0;
    case core::BondOrder::TRIPLE:
        return 3.0;
    case core::BondOrder::AROMATIC:
        return 1.5;
    }

    return 1.0;
}

[[nodiscard]] auto chi0(const core::Element& element) -> double {
    if (element.electronegativity <= 0.0) {
        throw std::logic_error{"MGC requires positive electronegativity for element '" +
                               std::string{element.symbol} + "'"};
    }

    return element.electronegativity;
}

} // namespace

auto MGCMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& table = core::periodic_table();

    const auto n = static_cast<Eigen::Index>(molecule.atom_count());

    if (n == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd X0 = Eigen::VectorXd::Zero(n);

    auto log_sum = 0.0;

    for (Eigen::Index i = 0; i < n; ++i) {
        const auto& atom = molecule.atom(static_cast<std::size_t>(i));
        const auto& element = table.element(atom.atomic_number());

        X0(i) = chi0(element);
        S(i, i) = 1.0;
        log_sum += std::log(X0(i));
    }

    for (const auto& bond : molecule.bonds()) {
        const auto i = static_cast<Eigen::Index>(bond.first_atom_index());
        const auto j = static_cast<Eigen::Index>(bond.second_atom_index());
        const auto w = order(bond.order());

        S(i, i) += w;
        S(j, j) += w;
        S(i, j) -= w;
        S(j, i) -= w;
    }

    Eigen::VectorXd chi = S.partialPivLu().solve(X0);
    chi -= X0;
    chi /= std::exp(log_sum / static_cast<double>(n));

    return charges::AtomicCharges{std::vector<double>{chi.data(), chi.data() + chi.size()}};
}

} // namespace chargefw::methods::builtin
