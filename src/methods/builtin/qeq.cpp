#include "methods/builtin/qeq.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/LU>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

using OverlapTerm = double (*)(double distance, double hardness_i, double hardness_j);

[[nodiscard]] auto resolve_overlap_term(const std::string& type) -> OverlapTerm {
    if (type == "Nishimoto-Mataga") {
        return [](const double distance, const double hardness_i, const double hardness_j) {
            return 1.0 / (distance + 2.0 / (hardness_i + hardness_j));
        };
    }
    if (type == "Nishimoto-Mataga-Weiss") {
        return [](const double distance, const double hardness_i, const double hardness_j) {
            constexpr auto factor = 1.2;
            return factor / (distance + (2.0 * factor) / (hardness_i + hardness_j));
        };
    }
    if (type == "Ohno") {
        return [](const double distance, const double hardness_i, const double hardness_j) {
            return 1.0 /
                   std::sqrt(distance * distance + std::pow(2.0 / (hardness_i + hardness_j), 2.0));
        };
    }
    if (type == "Ohno-Klopman") {
        return [](const double distance, const double hardness_i, const double hardness_j) {
            return 1.0 /
                   std::sqrt(distance * distance +
                             std::pow(1.0 / (2.0 * hardness_i) + 1.0 / (2.0 * hardness_j), 2.0));
        };
    }
    if (type == "DasGupta-Huzinaga") {
        return [](const double distance, const double hardness_i, const double hardness_j) {
            constexpr auto k = 0.4;
            return 1.0 / (distance + 1.0 / (hardness_i / 2.0 * std::exp(k * distance) +
                                            hardness_j / 2.0 * std::exp(k * distance)));
        };
    }
    if (type == "Louwen-Vogt") {
        return [](const double distance, const double hardness_i, const double hardness_j) {
            const auto gamma = (hardness_i + hardness_j) / 2.0;
            return 1.0 / std::cbrt(1.0 / std::pow(gamma, 3.0) + std::pow(distance, 3.0));
        };
    }

    throw std::logic_error{"Unsupported QEq overlap term '" + type + "'"};
}

} // namespace

auto QEqMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();
    const auto& geometry = input.geometry();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto electronegativity = parameters.atom("electronegativity");
    const auto hardness = parameters.atom("hardness");
    const auto& type = input.method_options().get<std::string>("overlap_term");

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        if (hardness[atom_index] <= 0.0) {
            throw std::logic_error{"QEq atom parameter 'hardness' must be positive"};
        }
    }

    const auto overlap_term = resolve_overlap_term(type);

    const auto n = static_cast<Eigen::Index>(atom_count);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);

        A(i, i) = hardness[atom_index];
        b(i) = -electronegativity[atom_index];

        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto value = overlap_term(geometry.distance(atom_index, other_atom_index),
                                            hardness[atom_index], hardness[other_atom_index]);

            A(i, j) = value;
            A(j, i) = value;
        }
    }

    A.row(n) = Eigen::VectorXd::Constant(n + 1, 1.0);
    A.col(n) = Eigen::VectorXd::Constant(n + 1, 1.0);
    A(n, n) = 0.0;
    b(n) = input.target_charge();

    const Eigen::VectorXd q = A.partialPivLu().solve(b).head(n);

    return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
}

} // namespace chargefw::methods::builtin
