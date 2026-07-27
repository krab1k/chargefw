#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cassert>
#include <cmath>
#include <vector>

namespace parameters = chargefw::parameters;

namespace {

auto make_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-kcm", .method_id = "kcm", .name = "Test KCM parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 1.0},
                                                    {.name = "hardness", .value = 1.0}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 2.0},
                                                    {.name = "hardness", .value = 1.0}}}}}};
}

auto assert_close(const double actual, const double expected) -> void {
    assert(std::abs(actual - expected) < 1.0e-12);
}

} // namespace

auto main() -> int {
    const auto charge_set = chargefw::test::calculate_single_method(chargefw::test::make_water(),
                                                                    "kcm", {make_parameter_set()});
    const auto& charges = charge_set.assignment(0).charges;

    chargefw::test::assert_calculation_provenance(charge_set, "kcm", "test-kcm");

    assert(charges.size() == 3);
    assert_close(charges[0], -0.4);
    assert_close(charges[1], 0.2);
    assert_close(charges[2], 0.2);
    assert_close(charges.total(), 0.0);

    return 0;
}
