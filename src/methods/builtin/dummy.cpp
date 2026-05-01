#include "methods/builtin/dummy.h"
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/methods/calculation_input.h>

namespace chargefw::methods::builtin {
auto DummyMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    return charges::AtomicCharges{std::vector(input.molecule().atom_count(), 0.0)};
}
} // namespace chargefw::methods::builtin
