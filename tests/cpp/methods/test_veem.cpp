#include "support/test_calculation.h"
#include "support/test_molecules.h"

#include <chargefw/methods/method_registry.h>

#include <cassert>

namespace methods = chargefw::methods;

auto main() -> int {
    const auto& registry = methods::method_registry();
    const auto* veem = registry.find("veem");

    assert(veem != nullptr);
    assert(!veem->requires_parameters());

    const auto workflow_charge_set =
        chargefw::test::calculate_method(chargefw::test::make_two_conformer_water(), "veem");
    chargefw::test::assert_calculation_provenance(workflow_charge_set, "veem", std::nullopt);
    chargefw::test::assert_conformer_independent(workflow_charge_set);

    return 0;
}
