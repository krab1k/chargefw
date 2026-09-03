#include "methods/builtin/sqeqp.h"

#include "methods/builtin/sqe.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace chargefw::methods::builtin {

auto SQEqpMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{sqe_core::calculate(input)};
    }

    const auto initial_charge_parameter = parameters.atom("q0");

    const auto n = static_cast<Eigen::Index>(atom_count);
    Eigen::VectorXd initial_charges = Eigen::VectorXd::Zero(n);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        initial_charges(static_cast<Eigen::Index>(atom_index)) =
            initial_charge_parameter[atom_index];
    }

    initial_charges.array() -= (initial_charges.sum() - core::total_formal_charge(molecule)) /
                               static_cast<double>(atom_count);

    return charges::AtomicCharges{sqe_core::calculate(input, {initial_charges.data(), atom_count})};
}

} // namespace chargefw::methods::builtin
