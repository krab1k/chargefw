#include "support/test_molecules.h"

#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>

#include <cassert>
#include <string>
#include <string_view>
#include <vector>

namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

auto calculate(const methods::Method& method, const chargefw::core::Molecule& molecule)
    -> chargefw::charges::AtomicCharges
{
    const features::PreparedMolecule prepared_molecule{molecule};
    const auto method_options = methods::make_default_options(method.option_schema());

    const methods::CalculationInput input{
        .prepared_molecule = prepared_molecule,
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

    const auto method_names = registry.names();
    assert((method_names == std::vector<std::string>{"dummy", "formal"}));

    assert(dummy->id() == std::string_view{"dummy"});
    assert(formal->id() == std::string_view{"formal"});

    assert(dummy->metadata().name == std::string_view{"Dummy method"});
    assert(dummy->metadata().full_name == std::string_view{"Dummy zero charges"});
    assert(!dummy->metadata().publication.has_value());
    assert(dummy->metadata().priority == 10);
    assert(!dummy->requirements().formal_charges);
    assert(dummy->option_schema().empty());

    assert(formal->metadata().name == std::string_view{"Formal"});
    assert(formal->metadata().full_name == std::string_view{"Formal atomic charges"});
    assert(!formal->metadata().publication.has_value());
    assert(formal->metadata().priority == 10);
    assert(formal->requirements().formal_charges);
    assert(formal->option_schema().empty());

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