#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
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
            .id = "test-eem", .method_id = "eem", .name = "Test EEM parameters"},
        parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
        parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 10.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 10.0}}}}}};
}

} // namespace

TEST_CASE("EEM produces conformer-dependent water charges", "[methods][eem]") {
    const auto charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "eem", {make_parameter_set()});
    const auto& charges = charge_set.assignment(0).charges;

    chargefw::test::assert_calculation_provenance(charge_set, "eem", "test-eem");
    chargefw::test::assert_conformer_dependent(charge_set, 2);

    CHECK(charges.size() == 3);
    CHECK(std::abs(charges[0] - (-0.07552338)) < 1.0e-8);
    CHECK(std::abs(charges[1] - (0.03776169)) < 1.0e-8);
    CHECK(std::abs(charges[2] - (0.03776169)) < 1.0e-8);
    CHECK(std::abs(charges.total() - (0.0)) < 1.0e-8);
    CHECK(std::abs(charges[0] - charge_set.assignment(1).charges[0]) > 1.0e-8);

    chargefw::test::assert_water_charges_labeling_invariant("eem", {make_parameter_set()});
}
