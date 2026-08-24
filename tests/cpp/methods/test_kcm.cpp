#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

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
            .id = "test-kcm", .method_id = "kcm", .name = "Test KCM parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 1.0},
                                                    {.name = "hardness", .value = 1.0}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 2.0},
                                                    {.name = "hardness", .value = 1.0}}}}}};
}

} // namespace

TEST_CASE("KCM produces conformer-independent water charges", "[methods][kcm]") {
    const auto charge_set = chargefw::test::calculate_single_method(
        chargefw::test::make_water_graph(), "kcm", {make_parameter_set()});
    const auto& charges = charge_set.assignment(0).charges;

    chargefw::test::assert_calculation_provenance(charge_set, "kcm", "test-kcm");
    chargefw::test::assert_conformer_independent(charge_set);

    CHECK(charges.size() == 3);
    CHECK(std::abs(charges[0] - (-0.4)) < 1.0e-12);
    CHECK(std::abs(charges[1] - (0.2)) < 1.0e-12);
    CHECK(std::abs(charges[2] - (0.2)) < 1.0e-12);
    CHECK(std::abs(charges.total() - (0.0)) < 1.0e-12);
}
