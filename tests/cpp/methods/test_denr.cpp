#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/methods/method_options.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <snitch/snitch.hpp>
#include <vector>

namespace parameters = chargefw::parameters;

namespace {

auto make_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-denr", .method_id = "denr", .name = "Test DENR parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 1.0},
                                                    {.name = "hardness", .value = 1.0}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 2.0},
                                                    {.name = "hardness", .value = 1.0}}}}}};
}

} // namespace

TEST_CASE("DENR produces conformer-independent water charges", "[methods][denr]") {
    auto options = chargefw::methods::MethodOptions{};
    options.set("step", 0.1);
    options.set("iterations", 100);

    const auto charge_set = chargefw::test::calculate_single_method(
        chargefw::test::make_water_graph(), "denr", {make_parameter_set()}, &options);
    const auto& charges = charge_set.assignment(0).charges;

    chargefw::test::assert_calculation_provenance(charge_set, "denr", "test-denr");
    chargefw::test::assert_conformer_independent(charge_set);
    chargefw::test::assert_neutral_water_charges(charges, 1.0e-12);

    chargefw::test::assert_water_charges_labeling_invariant("denr", {make_parameter_set()},
                                                            &options);
    chargefw::test::assert_water_charges_geometry_independent("denr", {make_parameter_set()},
                                                              &options);
}