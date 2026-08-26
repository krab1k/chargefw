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

TEST_CASE("EEM enforces a charged molecular target", "[methods][eem]") {
    const chargefw::core::Molecule cation{
        {chargefw::core::Atom{1, 1}},
        {},
        {chargefw::core::Conformer{{chargefw::core::Position{}}}}};
    const auto charge_set =
        chargefw::test::calculate_single_method(cation, "eem", {make_parameter_set()});

    CHECK(std::abs(charge_set.assignment(0).charges[0] - 1.0) < 1.0e-12);
    CHECK(std::abs(charge_set.assignment(0).charges.total() - 1.0) < 1.0e-12);
}

TEST_CASE("EEM has a stable water regression and responds to geometry", "[methods][eem]") {
    const auto charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "eem", {make_parameter_set()});
    const auto& charges = charge_set.assignment(0).charges;

    CHECK(std::abs(charges[0] - (-0.07552338)) < 1.0e-8);
    CHECK(std::abs(charges[1] - (0.03776169)) < 1.0e-8);
    CHECK(std::abs(charges[2] - (0.03776169)) < 1.0e-8);
    CHECK(std::abs(charges[0] - charge_set.assignment(1).charges[0]) > 1.0e-8);
}
