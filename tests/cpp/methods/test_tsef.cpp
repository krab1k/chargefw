#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>

#include <snitch/snitch.hpp>

namespace methods = chargefw::methods;

TEST_CASE("TSEF rejects unsupported unbonded molecules", "[methods][tsef]") {
    const auto* tsef = methods::method_registry().find("tsef");
    CHECK(tsef != nullptr);

    const chargefw::core::Molecule unbonded_pair{
        {chargefw::core::Atom{1}, chargefw::core::Atom{8}}};
    const chargefw::features::PreparedMolecule prepared_pair{unbonded_pair};
    const methods::MethodOptions options;
    const auto prerequisite_result = tsef->check_method_prerequisites(
        {.prepared_molecule = prepared_pair, .method_options = options});

    CHECK(!prerequisite_result);
    CHECK(prerequisite_result.issues().size() == 1);
    CHECK(prerequisite_result.issues()[0].kind ==
          methods::PrerequisiteIssueKind::unsupported_molecule);
}