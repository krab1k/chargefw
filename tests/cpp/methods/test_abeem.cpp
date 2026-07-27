#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cassert>
#include <cmath>
#include <string_view>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
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

auto assert_close(const double actual, const double expected) -> void {
    assert(std::abs(actual - expected) < 1.0e-4);
}

} // namespace

auto main() -> int {
    const auto collection =
        core::MoleculeCollection{std::vector{chargefw::test::make_water()}, "test"};
    const features::PreparedMoleculeCollection prepared{collection};

    const auto& registry = methods::method_registry();
    const auto* abeem = registry.find("abeem");

    assert(abeem != nullptr);

    const std::vector candidate_methods{abeem};
    const std::vector parameter_sets{make_parameter_set()};

    const auto applicability =
        methods::find_applicable_methods(prepared, candidate_methods, parameter_sets);

    assert(applicability.applicable.size() == 1);
    assert(applicability.rejected.empty());

    const auto charge_set = methods::calculate_charges(applicability.applicable.front(), prepared);
    const auto& charges = charge_set.assignment(0).charges;

    assert(charge_set.method_id() == std::string_view{"abeem"});
    assert(charge_set.parameter_set_id().has_value());
    assert(*charge_set.parameter_set_id() == std::string_view{"test-abeem"});

    assert(charges.size() == 3);
    assert(charges[0] < 0.0);
    assert(charges[1] > 0.0);
    assert(charges[2] > 0.0);
    assert_close(charges[1], charges[2]);
    assert_close(charges.total(), 0.0);

    return 0;
}
