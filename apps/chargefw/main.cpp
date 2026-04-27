#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/topology_features.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/methods/method_options.h>

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace chg = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

auto make_water() -> core::Molecule
{
    std::vector<core::Atom> atoms{
        core::Atom{8, 0, "O"},
        core::Atom{1, 0, "H1"},
        core::Atom{1, 0, "H2"}
    };

    std::vector<core::Bond> bonds{
        core::Bond{0, 1, core::BondOrder::SINGLE},
        core::Bond{0, 2, core::BondOrder::SINGLE}
    };

    std::vector<core::Position> positions{
        core::Position{0.0000, 0.0000, 0.0000},
        core::Position{0.9572, 0.0000, 0.0000},
        core::Position{-0.2390, 0.9270, 0.0000}
    };

    std::vector<core::Conformer> conformers{
        core::Conformer{std::move(positions), "model-1"}
    };

    return core::Molecule{
        std::move(atoms),
        std::move(bonds),
        std::move(conformers),
        "water"
    };
}

auto calculate_charges(const methods::Method& method, const core::Molecule& molecule)
    -> chg::AtomicCharges
{
    const features::TopologyFeatures topology{molecule};
    const features::ConformerFeatures geometry{molecule, 0};

    const methods::MethodOptions method_options{
        methods::make_default_options(method.option_schema())
    };

    const methods::CalculationInput input{
        .molecule = molecule,
        .topology = topology,
        .geometry = &geometry,
        .method_options = method_options
    };

    return method.calculate(input);
}

auto print_charges(std::string_view label,
                   const core::Molecule& molecule,
                   const chg::AtomicCharges& charges) -> void
{
    std::cout << label << '\n';

    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto& atom = molecule.atom(atom_index);

        std::cout
            << "  atom " << atom_index
            << " " << atom.name()
            << " Z=" << atom.atomic_number()
            << " q=" << charges[atom_index]
            << '\n';
    }

    std::cout << "  total = " << charges.total() << "\n\n";
}

auto require_method(const methods::MethodRegistry& registry, std::string_view id)
    -> const methods::Method&
{
    const auto* method = registry.find(id);

    if (method == nullptr) {
        throw std::runtime_error{"method not found: " + std::string{id}};
    }

    return *method;
}

auto main() -> int
{
    const auto molecule = make_water();

    const auto& registry = methods::method_registry();

    const auto& dummy = require_method(registry, "dummy");
    const auto& formal = require_method(registry, "formal");

    const auto dummy_charges = calculate_charges(dummy, molecule);
    const auto formal_charges = calculate_charges(formal, molecule);

    std::cout << "Molecule: " << molecule.name() << '\n';
    std::cout << "Atoms: " << molecule.atom_count() << '\n';
    std::cout << "Bonds: " << molecule.bond_count() << "\n\n";

    print_charges("Dummy charges:", molecule, dummy_charges);
    print_charges("Formal charges:", molecule, formal_charges);

    return 0;
}