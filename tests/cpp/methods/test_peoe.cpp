#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <chargefw/methods/method_options.h>
#include <cmath>
#include <snitch/snitch.hpp>

namespace parameters = chargefw::parameters;

namespace {

auto make_peoe_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "peoe-test-parameters", .method_id = "peoe", .name = "PEOE test parameters"},
        parameters::CommonParameters{{{.name = "dampH", .value = 20.02}}},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "A", .value = 7.17},
                                                    {.name = "B", .value = 6.24},
                                                    {.name = "C", .value = -0.56}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "A", .value = 12.06},
                                                    {.name = "B", .value = 13.85},
                                                    {.name = "C", .value = 3.98}}}}}};
}

} // namespace

TEST_CASE("PEOE iteration option changes water charges", "[methods][peoe]") {
    auto one_iteration = chargefw::methods::MethodOptions{};
    one_iteration.set("iters", 1);

    const auto default_charges = chargefw::test::calculate_single_method(
        chargefw::test::make_water_graph(), "peoe", {make_peoe_parameters()});
    const auto one_iteration_charges = chargefw::test::calculate_single_method(
        chargefw::test::make_water_graph(), "peoe", {make_peoe_parameters()}, &one_iteration);

    CHECK(std::abs(default_charges.assignment(0).charges[0] -
                   one_iteration_charges.assignment(0).charges[0]) > 1.0e-6);
}
