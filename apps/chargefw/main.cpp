#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/topology_features.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>

#include <cstddef>
#include <iostream>
#include <string_view>
#include <utility>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

auto print_molecule_summary(const core::Molecule& molecule,
                            const features::TopologyFeatures& topology,
                            const features::ConformerFeatures& geometry) -> void {
    std::cout << "Molecule: " << molecule.name() << '\n';
    std::cout << "Atoms: " << molecule.atom_count() << '\n';
    std::cout << "Bonds: " << molecule.bond_count() << "\n\n";

    std::cout << "Topology\n";
    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        std::cout << "  atom " << atom_index << " degree=" << topology.degree(atom_index)
                  << " neighbors=";

        for (const auto neighbor_index : topology.neighbor_indices(atom_index)) {
            std::cout << neighbor_index << ' ';
        }

        std::cout << '\n';
    }

    std::cout << "\nGeometry\n";
    std::cout << "  O-H1 distance = " << geometry.distance(0, 1) << '\n';
    std::cout << "  O-H2 distance = " << geometry.distance(0, 2) << "\n\n";
}

auto print_charges(const std::string_view method_id, const core::Molecule& molecule,
                   const chargefw::charges::AtomicCharges& charges) -> void {
    std::cout << "Method: " << method_id << '\n';

    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto& atom = molecule.atom(atom_index);

        std::cout << "  atom " << atom_index << ' ' << atom.name() << " Z=" << atom.atomic_number()
                  << " formal=" << atom.formal_charge() << " q=" << charges[atom_index] << '\n';
    }

    std::cout << "  total = " << charges.total() << "\n\n";
}

} // namespace

auto main() -> int {
    auto atoms = std::vector{
        core::Atom{8, 0, "O"},
        core::Atom{1, 0, "H1"},
        core::Atom{1, 0, "H2"},
    };

    auto bonds = std::vector{
        core::Bond{0, 1, core::BondOrder::SINGLE},
        core::Bond{0, 2, core::BondOrder::SINGLE},
    };

    auto positions = std::vector{
        core::Position{0.0000, 0.0000, 0.0000},
        core::Position{0.9572, 0.0000, 0.0000},
        core::Position{-0.2390, 0.9270, 0.0000},
    };

    auto conformers = std::vector{
        core::Conformer{std::move(positions), "model-1"},
    };

    const auto molecule = core::Molecule{
        std::move(atoms),
        std::move(bonds),
        std::move(conformers),
        "water",
    };

    const features::TopologyFeatures topology{molecule};
    const features::ConformerFeatures geometry{molecule, 0};

    print_molecule_summary(molecule, topology, geometry);

    const auto& registry = methods::method_registry();

    for (const auto& method_ptr : registry.methods()) {
        const auto& method = *method_ptr;
        const auto options = methods::make_default_options(method.option_schema());

        methods::validate_method_options(method.option_schema(), options);

        const auto input = methods::CalculationInput{
            .molecule = molecule,
            .topology = topology,
            .geometry = &geometry,
            .method_options = options,
        };

        std::cout << "Running method: " << method.id() << '\n';
        const auto charges = method.calculate(input);
        print_charges(method.id(), molecule, charges);
    }

    return 0;
}
