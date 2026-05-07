#include "support/test_molecules.h"

#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/models/atom_parameters.h>
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

auto atom_key(const int atomic_number) -> parameters::AtomParameterKey {
    return {.atomic_number = atomic_number,
            .classification = parameters::AtomParameterClassificationKind::PLAIN,
            .type = "*"};
}

auto make_parameter_set() -> parameters::ParameterSet {
    const parameters::ParameterSetMetadata metadata{
        .id = "test-denr", .method_id = "denr", .name = "Test DENR parameters"};
    const parameters::CommonParameters common{};
    const parameters::AtomParameters atom{
        {{.key = atom_key(1),
          .parameters = {{.name = "electronegativity", .value = 1.0},
                         {.name = "hardness", .value = 1.0}}},
         {.key = atom_key(8),
          .parameters = {{.name = "electronegativity", .value = 2.0},
                         {.name = "hardness", .value = 1.0}}}}};

    return parameters::ParameterSet{metadata, common, atom};
}

} // namespace

auto main() -> int {
    const auto collection =
        core::MoleculeCollection{std::vector{chargefw::test::make_water()}, "test"};
    const features::PreparedMoleculeCollection prepared{collection};

    const auto& registry = methods::method_registry();
    const auto* denr = registry.find("denr");

    assert(denr != nullptr);

    const std::vector<const methods::Method*> candidate_methods{denr};
    const std::vector parameter_sets{make_parameter_set()};

    auto applicability =
        methods::find_applicable_methods(prepared, candidate_methods, parameter_sets);

    assert(applicability.applicable.size() == 1);
    assert(applicability.rejected.empty());

    auto& selected = applicability.applicable.front();
    selected.method_options.set("step", 0.1);
    selected.method_options.set("iterations", 100);

    const auto charge_set = methods::calculate_charges(selected, prepared);
    const auto& charges = charge_set.assignment(0).charges;

    assert(charge_set.method_id() == std::string_view{"denr"});
    assert(charge_set.parameter_set_id().has_value());
    assert(*charge_set.parameter_set_id() == std::string_view{"test-denr"});

    assert(charges.size() == 3);
    assert(charges[0] < 0.0);
    assert(charges[1] > 0.0);
    assert(charges[2] > 0.0);
    assert(std::abs(charges[1] - charges[2]) < 1.0e-12);
    assert(std::abs(charges.total()) < 1.0e-12);

    return 0;
}