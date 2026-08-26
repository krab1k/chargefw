#include "support/test_calculation.h"
#include "support/test_molecules.h"

#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_registry.h>

#include <cmath>
#include <snitch/snitch.hpp>

namespace methods = chargefw::methods;

TEST_CASE("MGC rejects non-neutral molecules", "[methods][mgc]") {
    const auto& registry = methods::method_registry();
    const auto* mgc = registry.find("mgc");

    REQUIRE(mgc != nullptr);
    const chargefw::core::Molecule cation{{chargefw::core::Atom{1, 1}}};
    const chargefw::features::PreparedMolecule prepared_cation{cation};
    const auto prerequisite_result = mgc->check_method_prerequisites(
        {.prepared_molecule = prepared_cation, .method_options = {}});

    CHECK(!prerequisite_result);
    REQUIRE(prerequisite_result.issues().size() == 1);
    CHECK(prerequisite_result.issues()[0].kind ==
          methods::PrerequisiteIssueKind::unsupported_molecule);
}

TEST_CASE("MGC has a stable water graph-charge regression", "[methods][mgc]") {
    const auto charge_set =
        chargefw::test::calculate_single_method(chargefw::test::make_water_graph(), "mgc");
    const auto& charges = charge_set.assignment(0).charges;

    CHECK(std::abs(charges[0] - (-0.24280469939476498)) < 1.0e-12);
    CHECK(std::abs(charges[1] - 0.12140234969738249) < 1.0e-12);
    CHECK(std::abs(charges[2] - 0.12140234969738249) < 1.0e-12);
}
