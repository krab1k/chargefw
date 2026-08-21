#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cassert>
#include <cmath>
#include <vector>

namespace parameters = chargefw::parameters;

namespace {

auto make_mpeoe_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-mpeoe", .method_id = "mpeoe", .name = "Test MPEOE parameters"},
        parameters::CommonParameters{{{.name = "Hplus", .value = 3.435}}},
        parameters::AtomParameters{
            {{.key = chargefw::test::hbo_atom_key(1, "1"),
              .parameters = {{.name = "A", .value = 7.674}, {.name = "B", .value = 28.689}}},
             {.key = chargefw::test::hbo_atom_key(6, "1"),
              .parameters = {{.name = "A", .value = 7.976}, {.name = "B", .value = 4.625}}}}},
        parameters::BondParameters{{{.key = chargefw::test::plain_bond_key(),
                                     .parameters = {{.name = "f", .value = 0.505}}}}}};
}

} // namespace

auto main() -> int {
    const auto charge_set = chargefw::test::calculate_single_method(
        chargefw::test::make_methane_graph(), "mpeoe", {make_mpeoe_parameters()});

    chargefw::test::assert_calculation_provenance(charge_set, "mpeoe", "test-mpeoe");
    chargefw::test::assert_conformer_independent(charge_set);

    const auto& charges = charge_set.assignment(0).charges;

    assert(charges.size() == 5);

    const auto carbon_charge = charges[0];
    const auto hydrogen_charge = charges[1];

    assert(carbon_charge < 0.0);
    assert(hydrogen_charge > 0.0);

    assert(std::abs(charges[1] - charges[2]) < 1.0e-12);
    assert(std::abs(charges[1] - charges[3]) < 1.0e-12);
    assert(std::abs(charges[1] - charges[4]) < 1.0e-12);

    assert(std::abs(charges.total()) < 1.0e-12);

    return 0;
}
