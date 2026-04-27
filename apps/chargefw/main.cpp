#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/charges/charge_collection.h>


#include <iostream>
#include <utility>
#include <print>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;

auto make_water() -> core::Molecule
{
    std::vector<core::Atom> atoms{
        core::Atom{8, 0, "O"},
        core::Atom{1, 0 ,"H1"},
        core::Atom{1, 0,"H2"}
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

auto make_example_charges() -> chargefw::charges::ChargeCollection
{
    namespace chg = chargefw::charges;

    std::vector<chg::ChargeAssignment> peoe_assignments{
        chg::ChargeAssignment{
            .target = chg::ChargeTarget{
                .molecule_index = 0,
                .conformer_index = std::nullopt
            },
            .charges = chg::AtomicCharges{std::vector<double>{0.1, -0.2, 0.1}}
        }
    };

    std::vector<chg::ChargeAssignment> eem_assignments{
        chg::ChargeAssignment{
            .target = chg::ChargeTarget{
                .molecule_index = 0,
                .conformer_index = 0
            },
            .charges = chg::AtomicCharges{std::vector<double>{0.2, -0.3, 0.1}}
        }
    };

    return chg::ChargeCollection{
        std::vector<chg::ChargeSet>{
            chg::ChargeSet{
                "peoe",
                std::move(peoe_assignments)
            },
            chg::ChargeSet{
                "eem",
                std::move(eem_assignments),
                "parameters1"
            }
        }
    };
}

auto main() -> int
{

    auto collection = chargefw::core::MoleculeCollection{
        std::vector<chargefw::core::Molecule>{
            make_water(),
            make_water()
        },
        "example batch"
    };


    auto charges = make_example_charges();
    std::println("Charge set count: {}", charges.charge_set_count());
    std::println("{}", charges.charge_set(0).assignment(0).charges.values());


    const auto molecule = collection[0];
    const features::PreparedMolecule prepared{molecule};

    const features::ConformerFeatures geometry{molecule, 0};
    const auto oh_distance = geometry.distance(0, 1);
    const auto neighbors = geometry.neighbor_indices_within(0, 1.1);

    std::cout << oh_distance << std::endl;
    for (const auto &neighbor_index : prepared.neighbor_indices(1)) {
        std::cout << neighbor_index << ' ';
    }
    std::cout << '\n';

    std::cout << "Molecule: " << molecule.name() << '\n';
    std::cout << "Atoms: " << molecule.atom_count() << '\n';
    std::cout << "Bonds: " << molecule.bond_count() << '\n';
    std::cout << '\n';

    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto& atom = molecule.atom(atom_index);

        std::cout
            << "Atom " << atom_index
            << " name=" << atom.name()
            << " Z=" << atom.atomic_number()
            << " degree=" << prepared.degree(atom_index)
            << " neighbors=";

        for (const auto neighbor_index : prepared.neighbor_indices(atom_index)) {
            std::cout << neighbor_index << ' ';
        }

        std::cout << '\n';
    }

    std::cout << '\n';

    const auto oh_bond_index = prepared.bond_index_between(0, 1);

    if (oh_bond_index.has_value()) {
        const auto& bond = molecule.bond(*oh_bond_index);

        std::cout
            << "Atoms 0 and 1 are bonded. Bond index: "
            << *oh_bond_index
            << ", order: "
            << static_cast<int>(bond.order())
            << '\n';
    }

    std::cout
        << "Atoms 1 and 2 bonded? "
        << (prepared.are_bonded(1, 2) ? "yes" : "no")
        << '\n';



    return 0;
}