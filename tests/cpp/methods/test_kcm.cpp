#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <snitch/snitch.hpp>
#include <vector>

namespace parameters = chargefw::parameters;
namespace methods = chargefw::methods;

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

TEST_CASE("KCM rejects non-neutral molecules", "[methods][kcm]") {
    const auto& registry = methods::method_registry();
    const auto* kcm = registry.find("kcm");

    REQUIRE(kcm != nullptr);
    const chargefw::core::Molecule cation{{chargefw::core::Atom{1, 1}}};
    const chargefw::features::PreparedMolecule prepared_cation{cation};
    const auto prerequisite_result = kcm->check_method_prerequisites(
        {.prepared_molecule = prepared_cation, .method_options = {}});

    CHECK(!prerequisite_result);
    REQUIRE(prerequisite_result.issues().size() == 1);
    CHECK(prerequisite_result.issues()[0].kind ==
          methods::PrerequisiteIssueKind::unsupported_molecule);
}

TEST_CASE("KCM has a stable water regression", "[methods][kcm]") {
    const auto charge_set = chargefw::test::calculate_single_method(
        chargefw::test::make_water_graph(), "kcm", {make_parameter_set()});
    const auto& charges = charge_set.assignment(0).charges;

    CHECK(std::abs(charges[0] - (-0.4)) < 1.0e-12);
    CHECK(std::abs(charges[1] - (0.2)) < 1.0e-12);
    CHECK(std::abs(charges[2] - (0.2)) < 1.0e-12);
}
