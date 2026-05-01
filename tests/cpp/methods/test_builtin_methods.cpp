#include "support/test_molecules.h"

#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>

#include <cassert>
#include <cmath>
#include <string>
#include <string_view>
#include <vector>

namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

auto calculate(const methods::Method& method, const chargefw::core::Molecule& molecule)
    -> chargefw::charges::AtomicCharges {
    const features::PreparedMolecule prepared_molecule{molecule};
    const auto method_options = methods::make_default_options(method.option_schema());

    const methods::CalculationInput input{prepared_molecule, method_options};

    return method.calculate(input);
}

} // namespace

auto main() -> int {
    const auto& registry = methods::method_registry();

    const auto* dummy = registry.find("dummy");
    const auto* formal = registry.find("formal");
    const auto* veem = registry.find("veem");
    const auto* peoe = registry.find("peoe");

    assert(dummy != nullptr);
    assert(formal != nullptr);
    assert(veem != nullptr);
    assert(peoe != nullptr);

    const auto method_names = registry.names();
    assert((method_names == std::vector<std::string>{"dummy", "formal", "peoe", "veem"}));

    assert(dummy->id() == std::string_view{"dummy"});
    assert(formal->id() == std::string_view{"formal"});
    assert(veem->id() == std::string_view{"veem"});

    assert(dummy->metadata().name == std::string_view{"Dummy method"});
    assert(dummy->metadata().full_name == std::string_view{"Dummy zero charges"});
    assert(!dummy->metadata().publication.has_value());
    assert(dummy->metadata().priority == 10);
    assert(!dummy->requirements().formal_charges);
    assert(!dummy->requires_parameters());
    assert(dummy->option_schema().empty());

    assert(formal->metadata().name == std::string_view{"Formal"});
    assert(formal->metadata().full_name == std::string_view{"Formal atomic charges"});
    assert(!formal->metadata().publication.has_value());
    assert(formal->metadata().priority == 10);
    assert(formal->requirements().formal_charges);
    assert(!formal->requires_parameters());
    assert(formal->option_schema().empty());

    assert(veem->metadata().name == std::string_view{"VEEM"});
    assert(veem->metadata().full_name == std::string_view{"Valence Electrons Equalization Method"});
    assert(veem->metadata().publication.has_value());
    assert(veem->metadata().priority == 20);
    assert(veem->requirements().element_properties);
    assert(veem->requirements().resources.time == methods::ComplexityTerm::atoms);
    assert(veem->requirements().resources.memory == methods::ComplexityTerm::constant);
    assert(!veem->requires_parameters());
    assert(veem->option_schema().empty());

    assert(peoe->id() == std::string_view{"peoe"});
    assert(peoe->metadata().name == std::string_view{"PEOE"});
    assert(peoe->metadata().full_name ==
           std::string_view{"Partial Equalization of Atomic Electronegativity"});
    assert(peoe->metadata().publication.has_value());
    assert(peoe->metadata().priority == 120);

    assert(peoe->requirements().bond_graph);
    assert(peoe->requirements().requires_common_parameters());
    assert(peoe->requirements().requires_atom_parameters());
    assert(peoe->requirements().common_parameters.size() == 1);
    assert(peoe->requirements().atom_parameters.size() == 3);
    assert(peoe->requires_parameters());
    assert(peoe->option_schema().size() == 1);

    const auto water = chargefw::test::make_water();

    const auto dummy_charges = calculate(*dummy, water);
    assert(dummy_charges.size() == water.atom_count());

    for (const auto charge : dummy_charges.values()) {
        assert(charge == 0.0);
    }

    const auto veem_charges = calculate(*veem, water);
    assert(veem_charges.size() == water.atom_count());

    assert(veem_charges[0] < 0.0);
    assert(veem_charges[1] > 0.0);
    assert(veem_charges[2] > 0.0);
    assert(std::abs(veem_charges[1] - veem_charges[2]) < 1.0e-12);
    assert(std::abs(veem_charges.total()) < 1.0e-12);

    const auto charged_pair = chargefw::test::make_formally_charged_pair();

    const auto formal_charges = calculate(*formal, charged_pair);
    assert(formal_charges.size() == charged_pair.atom_count());

    assert(formal_charges[0] == 1.0);
    assert(formal_charges[1] == -1.0);
    assert(formal_charges.total() == 0.0);

    return 0;
}