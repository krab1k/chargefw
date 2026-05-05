#include "methods/builtin/eqeqc.h"

#include <chargefw/core/molecule.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/LU>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

constexpr auto lambda = 1.2;
constexpr auto k = 14.4;
constexpr auto hydrogen_atomic_number = 1;
constexpr auto hydrogen_electron_affinity = -2.0;

[[nodiscard]] auto electron_affinity_for_eqeq(const core::Element& element) -> double {
    if (element.atomic_number == hydrogen_atomic_number) {
        return hydrogen_electron_affinity;
    }

    return element.electron_affinity;
}

[[nodiscard]] auto interaction(const double distance, const double hardness_i,
                               const double hardness_j) -> double {

    if (hardness_i <= 0.0 || hardness_j <= 0.0) {
        throw std::logic_error{"EQeq+C hardness must be positive"};
    }

    const auto a = std::sqrt(hardness_i * hardness_j) / k;
    const auto overlap =
        std::exp(-a * a * distance * distance) * (2.0 * a - a * a * distance - 1.0 / distance);

    return lambda * k / 2.0 * (1.0 / distance + overlap);
}

[[nodiscard]] auto correction_term(const double alpha, const double distance,
                                   const double first_covalent_radius,
                                   const double second_covalent_radius, const double first_dz,
                                   const double second_dz) -> double {
    const auto tkk = first_dz - second_dz;
    const auto bkk = std::exp(-alpha * (distance - first_covalent_radius - second_covalent_radius));

    return tkk * bkk;
}

} // namespace

auto EQeqCMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& geometry = input.geometry();
    const auto& parameters = input.parameters();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto alpha = parameters.common("alpha");

    if (alpha <= 0.0) {
        throw std::logic_error{"EQeq+C common parameter 'alpha' must be positive"};
    }

    const auto dz = parameters.atom("Dz");

    const auto& table = core::periodic_table();
    const auto n = static_cast<Eigen::Index>(atom_count);

    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + 1);
    Eigen::VectorXd hardness = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd electronegativity = Eigen::VectorXd::Zero(n);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);
        const auto& atom = molecule.atom(atom_index);
        const auto& element = table.element(atom.atomic_number());
        const auto electron_affinity = electron_affinity_for_eqeq(element);

        electronegativity(i) = (element.first_ionization_potential + electron_affinity) / 2.0;
        hardness(i) = element.first_ionization_potential - electron_affinity;
    }

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);

        matrix(i, i) = hardness(i);
        rhs(i) = -electronegativity(i);

        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto distance = geometry.distance(atom_index, other_atom_index);
            const auto value = interaction(distance, hardness(i), hardness(j));

            matrix(i, j) = value;
            matrix(j, i) = value;
        }
    }

    matrix.row(n) = Eigen::VectorXd::Constant(n + 1, 1.0);
    matrix.col(n) = Eigen::VectorXd::Constant(n + 1, 1.0);
    matrix(n, n) = 0.0;
    rhs(n) = core::total_formal_charge(molecule);

    Eigen::VectorXd q = matrix.partialPivLu().solve(rhs).head(n);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto& atom = molecule.atom(atom_index);
        const auto& element = table.element(atom.atomic_number());

        auto correction = 0.0;

        for (std::size_t other_atom_index = 0; other_atom_index < atom_count; ++other_atom_index) {
            if (atom_index == other_atom_index) {
                continue;
            }

            const auto& other_atom = molecule.atom(other_atom_index);
            const auto& other_element = table.element(other_atom.atomic_number());
            const auto distance = geometry.distance(atom_index, other_atom_index);

            correction += correction_term(alpha, distance, element.covalent_radius,
                                          other_element.covalent_radius, dz[atom_index],
                                          dz[other_atom_index]);
        }

        q(static_cast<Eigen::Index>(atom_index)) += correction;
    }

    return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
}

} // namespace chargefw::methods::builtin