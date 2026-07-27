#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/methods/method_options.h>
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
    const parameters::ParameterSetMetadata metadata{
        .id = "test-denr", .method_id = "denr", .name = "Test DENR parameters"};
    const parameters::CommonParameters common{};
    const parameters::AtomParameters atom{
        {{.key = chargefw::test::plain_atom_key(1),
          .parameters = {{.name = "electronegativity", .value = 1.0},
                         {.name = "hardness", .value = 1.0}}},
         {.key = chargefw::test::plain_atom_key(8),
          .parameters = {{.name = "electronegativity", .value = 2.0},
                         {.name = "hardness", .value = 1.0}}}}};

    return parameters::ParameterSet{metadata, common, atom};
}

} // namespace

auto main() -> int {
    auto options = chargefw::methods::MethodOptions{};
    options.set("step", 0.1);
    options.set("iterations", 100);

    const auto charge_set = chargefw::test::calculate_single_method(
        chargefw::test::make_water(), "denr", {make_parameter_set()}, &options);
    const auto& charges = charge_set.assignment(0).charges;

    chargefw::test::assert_calculation_provenance(charge_set, "denr", "test-denr");

    assert(charges.size() == 3);
    assert(charges[0] < 0.0);
    assert(charges[1] > 0.0);
    assert(charges[2] > 0.0);
    assert(std::abs(charges[1] - charges[2]) < 1.0e-12);
    assert(std::abs(charges.total()) < 1.0e-12);

    return 0;
}
