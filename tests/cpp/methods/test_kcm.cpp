#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_registry.h>

#include <snitch/snitch.hpp>

namespace methods = chargefw::methods;

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
