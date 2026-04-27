#include "methods/builtin/formal.h"
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/methods/calculation_input.h>

namespace chargefw::methods::builtin {
[[nodiscard]] auto FormalMethod::calculate(const CalculationInput& input) const
    -> charges::AtomicCharges {

    std::vector<double> values;
    values.reserve(input.molecule.atom_count());

    for (const auto& atom : input.molecule.atoms()) {
        values.push_back(atom.formal_charge());
    }

    return charges::AtomicCharges{std::move(values)};
}

} // namespace chargefw::methods::builtin
