#include "support/test_calculation.h"
#include "support/test_molecules.h"

#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_registry.h>

#include <cmath>
#include <snitch/snitch.hpp>

namespace methods = chargefw::methods;

TEST_CASE("VEEM rejects unsupported elements", "[methods][veem]") {
    const auto& registry = methods::method_registry();
    const auto* veem = registry.find("veem");

    CHECK(veem != nullptr);
    const chargefw::core::Molecule neon_molecule{{chargefw::core::Atom{10}}};
    const chargefw::features::PreparedMolecule prepared_neon{neon_molecule};
    const auto prerequisite_result = veem->check_method_prerequisites(
        {.prepared_molecule = prepared_neon, .method_options = {}});

    CHECK(!prerequisite_result);
    CHECK(prerequisite_result.issues().size() == 1);
    CHECK(prerequisite_result.issues()[0].kind ==
          methods::PrerequisiteIssueKind::unsupported_molecule);
    CHECK(prerequisite_result.issues()[0].atom_index == 0);
}

TEST_CASE("VEEM has a stable water regression", "[methods][veem]") {
    const auto charge_set =
        chargefw::test::calculate_single_method(chargefw::test::make_water_graph(), "veem");
    const auto& charges = charge_set.assignment(0).charges;

    CHECK(std::abs(charges[0] - (-0.5942492012779552)) < 1.0e-12);
    CHECK(std::abs(charges[1] - 0.2971246006389776) < 1.0e-12);
    CHECK(std::abs(charges[2] - 0.2971246006389776) < 1.0e-12);
}
