#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>

#include <cassert>

namespace methods = chargefw::methods;

auto main() -> int {
    const auto* tsef = methods::method_registry().find("tsef");
    assert(tsef != nullptr);

    const chargefw::core::Molecule unbonded_pair{
        {chargefw::core::Atom{1}, chargefw::core::Atom{8}}};
    const chargefw::features::PreparedMolecule prepared_pair{unbonded_pair};
    const methods::MethodOptions options;
    const auto prerequisite_result = tsef->check_method_prerequisites(
        {.prepared_molecule = prepared_pair, .method_options = options});

    assert(!prerequisite_result);
    assert(prerequisite_result.issues().size() == 1);
    assert(prerequisite_result.issues()[0].kind ==
           methods::PrerequisiteIssueKind::unsupported_molecule);

    return 0;
}
