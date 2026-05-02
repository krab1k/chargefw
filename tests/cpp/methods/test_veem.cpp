#include "support/test_molecules.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/methods/method_prerequisites.h>

#include <cassert>
#include <stdexcept>
#include <vector>
#include <optional>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

auto make_iron_atom() -> core::Molecule {
    std::vector atoms{
        core::Atom{26, 0, "Fe"}
    };

    return core::Molecule{
        std::move(atoms),
        {},
        {},
        "iron"
    };
}

} // namespace

auto main() -> int {
    const auto& registry = methods::method_registry();
    const auto* veem = registry.find("veem");

    assert(veem != nullptr);
    assert(!veem->requires_parameters());

    const auto options = methods::make_default_options(veem->option_schema());

    const auto water = chargefw::test::make_water();
    const features::PreparedMolecule prepared_water{water};

    const auto water_check = veem->check_method_prerequisites({
        .prepared_molecule = prepared_water,
        .method_options = options
    });

    assert(water_check);

    const auto iron = make_iron_atom();
    const features::PreparedMolecule prepared_iron{iron};

    const auto iron_check = veem->check_method_prerequisites({
        .prepared_molecule = prepared_iron,
        .method_options = options
    });

    assert(!iron_check);
    assert(iron_check.issues().size() == 1);
    assert(iron_check.issues()[0].kind == methods::PrerequisiteIssueKind::unsupported_molecule);
    assert(iron_check.issues()[0].atom_index == std::optional<std::size_t>{0});

    bool rejected_unchecked_calculation = false;

    try {
        const methods::CalculationInput input{prepared_iron, options};

        [[maybe_unused]] const auto charges = veem->calculate(input);
    } catch (const std::logic_error&) {
        rejected_unchecked_calculation = true;
    }

    assert(rejected_unchecked_calculation);

    return 0;
}