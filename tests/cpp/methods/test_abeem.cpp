#include "support/test_assertions.h"
#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cassert>
#include <cmath>
#include <vector>

namespace parameters = chargefw::parameters;

namespace {

auto bond_key(const int first_atomic_number, const int second_atomic_number)
    -> parameters::BondParameterKey {
    return {.first_atom = chargefw::test::plain_atom_key(first_atomic_number),
            .second_atom = chargefw::test::plain_atom_key(second_atomic_number),
            .bond = {.classification = parameters::BondParameterClassificationKind::PLAIN,
                     .type = "*"}};
}

auto make_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-abeem", .method_id = "abeem", .name = "Test ABEEM parameters"},
        parameters::CommonParameters{{{.name = "k", .value = 1.0}}},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "a", .value = 1.0},
                                                    {.name = "b", .value = 10.0},
                                                    {.name = "c", .value = 0.5}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "a", .value = 2.0},
                                                    {.name = "b", .value = 10.0},
                                                    {.name = "c", .value = 0.5}}}}},
        parameters::BondParameters{{{.key = bond_key(8, 1),
                                     .parameters = {{.name = "A", .value = 1.0},
                                                    {.name = "B", .value = 10.0},
                                                    {.name = "C", .value = 0.5},
                                                    {.name = "D", .value = 0.5}}}}}};
}

} // namespace

auto main() -> int {
    const auto charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "abeem", {make_parameter_set()});
    const auto& charges = charge_set.assignment(0).charges;

    chargefw::test::assert_calculation_provenance(charge_set, "abeem", "test-abeem");
    chargefw::test::assert_conformer_dependent(charge_set, 2);

    assert(charges.size() == 3);
    assert(charges[0] < 0.0);
    assert(charges[1] > 0.0);
    assert(charges[2] > 0.0);
    chargefw::test::assert_close(charges[1], charges[2], 1.0e-4);
    chargefw::test::assert_close(charges.total(), 0.0, 1.0e-4);
    assert(std::abs(charges[0] - charge_set.assignment(1).charges[0]) > 1.0e-4);

    return 0;
}
