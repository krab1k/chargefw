#pragma once

#include <Eigen/Core>
#include <Eigen/LU>

namespace chargefw::methods::builtin {

// Solves matrix * x = rhs with partial-pivoting LU. The matrix is overwritten by its factors
// instead of copied, which halves peak memory and cache footprint for large dense systems.
[[nodiscard]] inline auto solve_in_place(Eigen::MatrixXd& matrix, const Eigen::VectorXd& rhs)
    -> Eigen::VectorXd {
    const Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXd>> lu{matrix};
    return lu.solve(rhs);
}

} // namespace chargefw::methods::builtin
