#include "support/test_assertions.h"
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
            .id = "test-qeq", .method_id = "qeq", .name = "Test QEq parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 4.5280},
                                                    {.name = "hardness", .value = 13.8904}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 8.741},
                                                    {.name = "hardness", .value = 13.364}}}}}};
}

} // namespace

auto main() -> int {
    const auto charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "qeq", {make_parameter_set()});
    const auto& charges = charge_set.assignment(0).charges;

    chargefw::test::assert_calculation_provenance(charge_set, "qeq", "test-qeq");
    chargefw::test::assert_conformer_dependent(charge_set, 2);

    assert(charges.size() == 3);
    chargefw::test::assert_close(charges[0], -0.22710866, 1.0e-8);
    chargefw::test::assert_close(charges[1], 0.11355540, 1.0e-8);
    chargefw::test::assert_close(charges[2], 0.11355326, 1.0e-8);
    chargefw::test::assert_close(charges.total(), 0.0, 1.0e-8);
    assert(std::abs(charges[0] - charge_set.assignment(1).charges[0]) > 1.0e-8);

    return 0;
}
