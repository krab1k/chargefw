#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cmath>
#include <snitch/snitch.hpp>
#include <vector>

namespace parameters = chargefw::parameters;

namespace {

auto make_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-sqeqp", .method_id = "sqeqp", .name = "Test SQE+qp parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 4.5280},
                                                    {.name = "hardness", .value = 13.8904},
                                                    {.name = "width", .value = 1.0},
                                                    {.name = "q0", .value = 0.25}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 8.741},
                                                    {.name = "hardness", .value = 13.364},
                                                    {.name = "width", .value = 1.0},
                                                    {.name = "q0", .value = -0.5}}}}},
        parameters::BondParameters{{{.key = chargefw::test::single_bond_key(1, 8),
                                     .parameters = {{.name = "kappa", .value = 1.0}}}}}};
}

} // namespace

TEST_CASE("SQE+qp produces conformer-dependent water charges", "[methods][sqeqp]") {
    const auto charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "sqeqp", {make_parameter_set()});
    const auto& charges = charge_set.assignment(0).charges;

    chargefw::test::assert_calculation_provenance(charge_set, "sqeqp", "test-sqeqp");
    chargefw::test::assert_conformer_dependent(charge_set, 2);

    chargefw::test::assert_neutral_water_charges(charges, 1.0e-4);
    CHECK(std::abs(charges.total() - (0.0)) < 1.0e-12);
    CHECK(std::abs(charges[0] - charge_set.assignment(1).charges[0]) > 1.0e-4);

    chargefw::test::assert_water_charges_labeling_invariant("sqeqp", {make_parameter_set()});
}