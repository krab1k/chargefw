#include "methods/builtin/tsef.h"

#include "features/topology_helpers.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/LU>

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

constexpr auto alpha = 14.4;

[[nodiscard]] auto k(const int distance) -> double {
    constexpr std::array values{0.556, 0.778, 1.000, 1.053, 1.087, 1.091};

    if (std::cmp_greater(distance, values.size())) {
        return values.back();
    }

    return values.at(static_cast<std::size_t>(distance - 1));
}

[[nodiscard]] auto interaction(const int bond_distance) -> double {
    return alpha * k(bond_distance) / (0.84 * static_cast<double>(bond_distance) + 0.46);
}

} // namespace

auto TSEFMethod::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                         PrerequisiteResult& result) const -> void {
    if (!features::is_connected(input.prepared_molecule.topology().adjacency())) {
        result.add(PrerequisiteIssue{.kind = PrerequisiteIssueKind::unsupported_molecule,
                                     .message = "TSEF requires a connected molecular bond graph"});
    }
}

auto TSEFMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto electronegativity = parameters.atom("electronegativity");
    const auto hardness = parameters.atom("hardness");

    const auto n = static_cast<Eigen::Index>(atom_count);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);

    const auto bond_distances = features::all_pairs_bond_distances(input.topology().adjacency());

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);

        A(i, i) = hardness[atom_index];
        b(i) = -electronegativity[atom_index];

        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto bond_distance = bond_distances[atom_index][other_atom_index];
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto value = interaction(bond_distance);

            A(i, j) = value;
            A(j, i) = value;
        }
    }

    A.row(n) = Eigen::VectorXd::Constant(n + 1, 1.0);
    A.col(n) = Eigen::VectorXd::Constant(n + 1, 1.0);
    A(n, n) = 0.0;
    b(n) = core::total_formal_charge(molecule);

    const Eigen::VectorXd q = A.partialPivLu().solve(b).head(n);

    return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
}

} // namespace chargefw::methods::builtin
