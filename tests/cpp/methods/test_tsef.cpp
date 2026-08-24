#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cmath>
#include <snitch/snitch.hpp>

namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

auto tsef_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-tsef", .method_id = "tsef", .name = "Test TSEF parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 4.5280},
                                                    {.name = "hardness", .value = 13.8904}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 8.741},
                                                    {.name = "hardness", .value = 13.364}}}}}};
}

} // namespace

TEST_CASE("TSEF rejects unsupported unbonded molecules", "[methods][tsef]") {
    const auto* tsef = methods::method_registry().find("tsef");
    REQUIRE(tsef != nullptr);

    const chargefw::core::Molecule unbonded_pair{
        {chargefw::core::Atom{1}, chargefw::core::Atom{8}}};
    const chargefw::features::PreparedMolecule prepared_pair{unbonded_pair};
    const methods::MethodOptions options;
    const auto prerequisite_result = tsef->check_method_prerequisites(
        {.prepared_molecule = prepared_pair, .method_options = options});

    CHECK(!prerequisite_result);
    REQUIRE(prerequisite_result.issues().size() == 1);
    CHECK(prerequisite_result.issues()[0].kind ==
          methods::PrerequisiteIssueKind::unsupported_molecule);
}

TEST_CASE("TSEF has a stable water regression", "[methods][tsef]") {
    const auto charge_set = chargefw::test::calculate_single_method(
        chargefw::test::make_water_graph(), "tsef", {tsef_parameters()});
    const auto& charges = charge_set.assignment(0).charges;

    CHECK(std::abs(charges[0] - (-0.3971069862707636)) < 1.0e-12);
    CHECK(std::abs(charges[1] - 0.1985534931353818) < 1.0e-12);
    CHECK(std::abs(charges[2] - 0.1985534931353818) < 1.0e-12);
}
