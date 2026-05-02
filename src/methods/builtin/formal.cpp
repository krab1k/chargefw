#include "methods/builtin/formal.h"
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/methods/calculation_input.h>

namespace chargefw::methods::builtin {
auto FormalMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    std::vector<double> values;

    const auto& molecule = input.molecule();
    values.reserve(molecule.atom_count());

    for (const auto& atom : molecule.atoms()) {
        values.push_back(atom.formal_charge());
    }

    return charges::AtomicCharges{std::move(values)};
}

} // namespace chargefw::methods::builtin
