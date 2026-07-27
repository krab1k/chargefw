#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
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
            .id = "test-eqeqc", .method_id = "eqeqc", .name = "Test EQeq+C parameters"},
        parameters::CommonParameters{{{.name = "alpha", .value = 1.0}}},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "Dz", .value = 0.1}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "Dz", .value = 0.2}}}}}};
}

auto assert_close(const double actual, const double expected) -> void {
    assert(std::abs(actual - expected) < 1.0e-8);
}

} // namespace

auto main() -> int {
    const auto charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "eqeqc", {make_parameter_set()});
    const auto& charges = charge_set.assignment(0).charges;

    chargefw::test::assert_calculation_provenance(charge_set, "eqeqc", "test-eqeqc");
    chargefw::test::assert_conformer_dependent(charge_set, 2);

    assert(charges.size() == 3);
    assert_close(charges[0], -0.16494533);
    assert_close(charges[1], 0.08248506);
    assert_close(charges[2], 0.08246027);
    assert_close(charges.total(), 0.0);
    assert(std::abs(charges[0] - charge_set.assignment(1).charges[0]) > 1.0e-8);

    return 0;
}
