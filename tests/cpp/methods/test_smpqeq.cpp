#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
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
            .id = "test-smpqeq", .method_id = "smpqeq", .name = "Test SMP/QEq parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "first", .value = 1.0},
                                                    {.name = "second", .value = 10.0},
                                                    {.name = "third", .value = 1.0},
                                                    {.name = "fourth", .value = 0.1}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "first", .value = 2.0},
                                                    {.name = "second", .value = 10.0},
                                                    {.name = "third", .value = 1.0},
                                                    {.name = "fourth", .value = 0.1}}}}}};
}

} // namespace

TEST_CASE("SMP/QEq responds to geometry", "[methods][smpqeq]") {
    const auto charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "smpqeq", {make_parameter_set()});

    CHECK(std::abs(charge_set.assignment(0).charges[0] - charge_set.assignment(1).charges[0]) >
          1.0e-8);
}
