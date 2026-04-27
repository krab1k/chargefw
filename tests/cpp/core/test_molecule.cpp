#include "support/test_molecules.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>

#include <cassert>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace core = chargefw::core;

auto main() -> int
{
    const auto water = chargefw::test::make_water();

    assert(water.name() == std::string_view{"water"});
    assert(water.atom_count() == 3);
    assert(water.bond_count() == 2);
    assert(water.conformer_count() == 1);
    assert(water.has_coordinates());

    bool rejected_invalid_bond_index = false;

    try {
        std::vector atoms{
            core::Atom{6, 0, "C"}
        };

        std::vector bonds{
            core::Bond{0, 1, core::BondOrder::SINGLE}
        };

        [[maybe_unused]] const core::Molecule invalid{
            std::move(atoms),
            std::move(bonds)
        };
    } catch (const std::invalid_argument&) {
        rejected_invalid_bond_index = true;
    }

    assert(rejected_invalid_bond_index);

    bool rejected_duplicate_bond = false;

    try {
        std::vector atoms{
            core::Atom{6, 0, "C"},
            core::Atom{1, 0, "H"}
        };

        std::vector bonds{
            core::Bond{0, 1, core::BondOrder::SINGLE},
            core::Bond{1, 0, core::BondOrder::SINGLE}
        };

        [[maybe_unused]] const core::Molecule invalid{
            std::move(atoms),
            std::move(bonds)
        };
    } catch (const std::invalid_argument&) {
        rejected_duplicate_bond = true;
    }

    assert(rejected_duplicate_bond);

    bool rejected_wrong_conformer_size = false;

    try {
        std::vector atoms{
            core::Atom{8, 0,"O"},
            core::Atom{1, 0, "H1"}
        };

        std::vector positions{
            core::Position{0.0, 0.0, 0.0}
        };

        std::vector conformers{
            core::Conformer{std::move(positions)}
        };

        [[maybe_unused]] const core::Molecule invalid{
            std::move(atoms),
            {},
            std::move(conformers)
        };
    } catch (const std::invalid_argument&) {
        rejected_wrong_conformer_size = true;
    }

    assert(rejected_wrong_conformer_size);

    return 0;
}