#include "support/test_calculation.h"
#include "support/test_molecules.h"

#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_registry.h>

#include <cassert>

namespace methods = chargefw::methods;

auto main() -> int {
    const auto& registry = methods::method_registry();
    const auto* veem = registry.find("veem");

    assert(veem != nullptr);
    assert(!veem->requires_parameters());

    const chargefw::core::Molecule neon_molecule{{chargefw::core::Atom{10}}};
    const chargefw::features::PreparedMolecule prepared_neon{neon_molecule};
    const auto prerequisite_result = veem->check_method_prerequisites(
        {.prepared_molecule = prepared_neon, .method_options = {}});

    assert(!prerequisite_result);
    assert(prerequisite_result.issues().size() == 1);
    assert(prerequisite_result.issues()[0].kind ==
           methods::PrerequisiteIssueKind::unsupported_molecule);
    assert(prerequisite_result.issues()[0].atom_index == 0);

    const auto workflow_charge_set =
        chargefw::test::calculate_method(chargefw::test::make_two_conformer_water(), "veem");
    chargefw::test::assert_calculation_provenance(workflow_charge_set, "veem", std::nullopt);
    chargefw::test::assert_conformer_independent(workflow_charge_set);

    chargefw::test::assert_water_charges_geometry_independent("veem");
    chargefw::test::assert_water_charges_labeling_invariant("veem");

    return 0;
}
