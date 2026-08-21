#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cassert>
#include <cmath>

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
                                    {.key = chargefw::test::plain_atom_key(9),
                                     .parameters = {{.name = "A", .value = 12.06},
                                                    {.name = "B", .value = 13.85},
                                                    {.name = "C", .value = 3.98}}}}}};
}

} // namespace

auto main() -> int {
    const auto charge_set = chargefw::test::calculate_single_method(
        chargefw::test::make_hf_graph(), "peoe", {make_peoe_parameters()});

    chargefw::test::assert_calculation_provenance(charge_set, "peoe", "peoe-test-parameters");
    chargefw::test::assert_conformer_independent(charge_set);

    const auto& charges = charge_set.assignment(0).charges;

    assert(charges.size() == 2);
    assert(charges[0] > 0.0);
    assert(charges[1] < 0.0);
    assert(std::abs(charges.total()) < 1.0e-12);

    return 0;
}
