#include "support/test_molecules.h"

#include <chargefw/features/topology_features.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>

#include <cassert>
#include <string_view>

namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

auto calculate(const methods::Method& method, const chargefw::core::Molecule& molecule)
    -> chargefw::charges::AtomicCharges
{
    const features::TopologyFeatures topology{molecule};
    const auto method_options = methods::make_default_options(method.option_schema());

    const methods::CalculationInput input{
        .molecule = molecule,
        .topology = topology,
        .geometry = nullptr,
        .method_options = method_options
    };

    return method.calculate(input);
}

} // namespace

auto main() -> int
{
    const auto& registry = methods::method_registry();

    const auto* dummy = registry.find("dummy");
    const auto* formal = registry.find("formal");

    assert(dummy != nullptr);
    assert(formal != nullptr);

    assert(dummy->id() == std::string_view{"dummy"});
    assert(formal->id() == std::string_view{"formal"});

    const auto water = chargefw::test::make_water();

    const auto dummy_charges = calculate(*dummy, water);
    assert(dummy_charges.size() == water.atom_count());

    for (const auto charge : dummy_charges.values()) {
        assert(charge == 0.0);
    }

    const auto charged_pair = chargefw::test::make_formally_charged_pair();

    const auto formal_charges = calculate(*formal, charged_pair);
    assert(formal_charges.size() == charged_pair.atom_count());

    assert(formal_charges[0] == 1.0);
    assert(formal_charges[1] == -1.0);
    assert(formal_charges.total() == 0.0);

    return 0;
}