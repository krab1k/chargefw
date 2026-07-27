#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/models/atom_parameters.h>
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

auto assert_close(const double actual, const double expected) -> void {
    assert(std::abs(actual - expected) < 1.0e-8);
}

} // namespace

auto main() -> int {
    const auto collection =
        core::MoleculeCollection{std::vector{chargefw::test::make_water()}, "test"};
    const features::PreparedMoleculeCollection prepared{collection};

    const auto& registry = methods::method_registry();
    const auto* smpqeq = registry.find("smpqeq");

    assert(smpqeq != nullptr);

    const std::vector candidate_methods{smpqeq};
    const std::vector parameter_sets{make_parameter_set()};

    const auto applicability =
        methods::find_applicable_methods(prepared, candidate_methods, parameter_sets);

    assert(applicability.applicable.size() == 1);
    assert(applicability.rejected.empty());

    const auto charge_set = methods::calculate_charges(applicability.applicable.front(), prepared);
    const auto& charges = charge_set.assignment(0).charges;

    assert(charge_set.method_id() == std::string_view{"smpqeq"});
    assert(charge_set.parameter_set_id().has_value());
    assert(*charge_set.parameter_set_id() == std::string_view{"test-smpqeq"});

    assert(charges.size() == 3);
    assert_close(charges[0], -0.035475768417);
    assert_close(charges[1], 0.017737997826);
    assert_close(charges[2], 0.017737770591);
    assert_close(charges.total(), 0.0);

    return 0;
}
