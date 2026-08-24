#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cmath>
#include <snitch/snitch.hpp>
#include <vector>

namespace parameters = chargefw::parameters;

namespace {

auto make_mpeoe_parameters(const double bond_attenuation) -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-mpeoe", .method_id = "mpeoe", .name = "Test MPEOE parameters"},
        parameters::CommonParameters{{{.name = "Hplus", .value = 3.435}}},
        parameters::AtomParameters{
            {{.key = chargefw::test::hbo_atom_key(1, "1"),
              .parameters = {{.name = "A", .value = 7.674}, {.name = "B", .value = 28.689}}},
             {.key = chargefw::test::hbo_atom_key(8, "1"),
              .parameters = {{.name = "A", .value = 7.976}, {.name = "B", .value = 4.625}}}}},
        parameters::BondParameters{{{.key = chargefw::test::plain_bond_key(),
                                     .parameters = {{.name = "f", .value = bond_attenuation}}}}}};
}

} // namespace

TEST_CASE("MPEOE bond attenuation changes water charges", "[methods][mpeoe]") {
    const auto reference_charges = chargefw::test::calculate_single_method(
        chargefw::test::make_water_graph(), "mpeoe", {make_mpeoe_parameters(0.505)});
    const auto lower_attenuation_charges = chargefw::test::calculate_single_method(
        chargefw::test::make_water_graph(), "mpeoe", {make_mpeoe_parameters(0.25)});

    CHECK(std::abs(reference_charges.assignment(0).charges[0] -
                   lower_attenuation_charges.assignment(0).charges[0]) > 1.0e-6);
}
