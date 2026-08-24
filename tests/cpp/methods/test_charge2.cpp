#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/methods/method_options.h>
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

auto charge2_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-charge2", .method_id = "charge2", .name = "Test Charge2 parameters"},
        parameters::CommonParameters{{{.name = "a1", .value = 40.4},
                                      {.name = "a2", .value = 23.8},
                                      {.name = "a3", .value = 21.0},
                                      {.name = "b", .value = 198.4},
                                      {.name = "c", .value = 5.0},
                                      {.name = "alpha", .value = 3.0}}},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "chi", .value = 7.17},
                                                    {.name = "P0", .value = 0.628},
                                                    {.name = "q0", .value = 0.034}}},
                                    {.key = chargefw::test::plain_atom_key(6),
                                     .parameters = {{.name = "chi", .value = 7.98},
                                                    {.name = "P0", .value = 1.056},
                                                    {.name = "q0", .value = -0.068}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "chi", .value = 15.25},
                                                    {.name = "P0", .value = 0.912},
                                                    {.name = "q0", .value = -0.350}}}}}};
}

auto make_hco_chain() -> chargefw::core::Molecule {
    return chargefw::core::Molecule{{chargefw::core::Atom{1, 0, "H"},
                                     chargefw::core::Atom{6, 0, "C"},
                                     chargefw::core::Atom{8, 0, "O"}},
                                    {chargefw::core::Bond{0, 1}, chargefw::core::Bond{1, 2}},
                                    {},
                                    "HCO chain"};
}

} // namespace

TEST_CASE("Charge2 has a stable water electronegativity-transfer regression",
          "[methods][charge2]") {
    const auto charge_set = chargefw::test::calculate_single_method(
        chargefw::test::make_water_graph(), "charge2", {charge2_parameters()});
    const auto& charges = charge_set.assignment(0).charges;

    CHECK(std::abs(charges[0] - (-0.6789915966386555)) < 1.0e-12);
    CHECK(std::abs(charges[1] - 0.3394957983193277) < 1.0e-12);
    CHECK(std::abs(charges[2] - 0.3394957983193277) < 1.0e-12);
}

TEST_CASE("Charge2 iteration option changes polarizability feedback", "[methods][charge2]") {
    auto one_iteration = chargefw::methods::MethodOptions{};
    one_iteration.set("iters", 1);

    const auto default_charges = chargefw::test::calculate_single_method(
        make_hco_chain(), "charge2", {charge2_parameters()});
    const auto one_iteration_charges = chargefw::test::calculate_single_method(
        make_hco_chain(), "charge2", {charge2_parameters()}, &one_iteration);

    CHECK(std::abs(default_charges.assignment(0).charges[0] -
                   one_iteration_charges.assignment(0).charges[0]) > 1.0e-6);
}
