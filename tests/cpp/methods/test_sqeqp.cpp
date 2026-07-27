#include "support/test_molecules.h"

#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
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

auto atom_key(const int atomic_number) -> parameters::AtomParameterKey {
    return {.atomic_number = atomic_number,
            .classification = parameters::AtomParameterClassificationKind::PLAIN,
            .type = "*"};
}

auto bond_key(const int first_atomic_number, const int second_atomic_number)
    -> parameters::BondParameterKey {
    return {.first_atom = atom_key(first_atomic_number),
            .second_atom = atom_key(second_atomic_number),
            .bond = {.classification = parameters::BondParameterClassificationKind::BOND_ORDER,
                     .type = "1"}};
}

auto make_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-sqeqp", .method_id = "sqeqp", .name = "Test SQE+qp parameters"},
        {},
        parameters::AtomParameters{{{.key = atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 4.5280},
                                                    {.name = "hardness", .value = 13.8904},
                                                    {.name = "width", .value = 1.0},
                                                    {.name = "q0", .value = 0.25}}},
                                    {.key = atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 8.741},
                                                    {.name = "hardness", .value = 13.364},
                                                    {.name = "width", .value = 1.0},
                                                    {.name = "q0", .value = -0.5}}}}},
        parameters::BondParameters{
            {{.key = bond_key(1, 8), .parameters = {{.name = "kappa", .value = 1.0}}}}}};
}

} // namespace

auto main() -> int {
    const auto collection =
        core::MoleculeCollection{std::vector{chargefw::test::make_water()}, "test"};
    const features::PreparedMoleculeCollection prepared{collection};

    const auto& registry = methods::method_registry();
    const auto* sqeqp = registry.find("sqeqp");

    assert(sqeqp != nullptr);

    const std::vector candidate_methods{sqeqp};
    const std::vector parameter_sets{make_parameter_set()};

    const auto applicability =
        methods::find_applicable_methods(prepared, candidate_methods, parameter_sets);

    assert(applicability.applicable.size() == 1);
    assert(applicability.rejected.empty());

    const auto charge_set = methods::calculate_charges(applicability.applicable.front(), prepared);
    const auto& charges = charge_set.assignment(0).charges;

    assert(charge_set.method_id() == std::string_view{"sqeqp"});
    assert(charge_set.parameter_set_id().has_value());
    assert(*charge_set.parameter_set_id() == std::string_view{"test-sqeqp"});

    assert(charges.size() == 3);
    assert(charges[0] < 0.0);
    assert(charges[1] > 0.0);
    assert(charges[2] > 0.0);
    assert(std::abs(charges[1] - charges[2]) < 1.0e-4);
    assert(std::abs(charges.total()) < 1.0e-12);

    return 0;
}
