#include "support/test_molecules.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/position.h>

#include <string>
#include <utility>
#include <vector>

namespace chargefw::test {
namespace {

auto make_water_graph_parts() -> std::pair<std::vector<core::Atom>, std::vector<core::Bond>> {
    std::vector atoms{core::Atom{8, 0, "O"}, core::Atom{1, 0, "H1"}, core::Atom{1, 0, "H2"}};

    std::vector bonds{core::Bond{0, 1, core::BondOrder::SINGLE},
                      core::Bond{0, 2, core::BondOrder::SINGLE}};

    return {std::move(atoms), std::move(bonds)};
}

auto make_water_positions() -> std::vector<core::Position> {
    return {core::Position{.x = 0.0000, .y = 0.0000, .z = 0.0000},
            core::Position{.x = 0.7570, .y = 0.5859, .z = 0.0000},
            core::Position{.x = -0.7570, .y = 0.5859, .z = 0.0000}};
}

} // namespace

auto make_water_graph() -> core::Molecule {
    auto [atoms, bonds] = make_water_graph_parts();

    return core::Molecule{std::move(atoms), std::move(bonds), {}, "water"};
}

auto make_water() -> core::Molecule {
    auto [atoms, bonds] = make_water_graph_parts();

    std::vector conformers{core::Conformer{make_water_positions(), "model-1"}};

    return core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers), "water"};
}

auto make_two_conformer_water() -> core::Molecule {
    auto [atoms, bonds] = make_water_graph_parts();

    std::vector positions_2{core::Position{.x = 0.0000, .y = 0.0000, .z = 0.0000},
                            core::Position{.x = 1.1000, .y = 0.0000, .z = 0.0000},
                            core::Position{.x = -0.3000, .y = 1.0500, .z = 0.0000}};

    std::vector conformers{core::Conformer{make_water_positions(), "model-1"},
                           core::Conformer{std::move(positions_2), "model-2"}};

    return core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers), "water"};
}

auto make_hf_graph() -> core::Molecule {
    std::vector atoms{core::Atom{1, 0, "H"}, core::Atom{9, 0, "F"}};

    std::vector bonds{core::Bond{0, 1, core::BondOrder::SINGLE}};

    return core::Molecule{std::move(atoms), std::move(bonds), {}, "hf"};
}

auto make_methane_graph() -> core::Molecule {
    std::vector atoms{core::Atom{6, 0, "C"}, core::Atom{1, 0, "H1"}, core::Atom{1, 0, "H2"},
                      core::Atom{1, 0, "H3"}, core::Atom{1, 0, "H4"}};

    std::vector bonds{
        core::Bond{0, 1, core::BondOrder::SINGLE}, core::Bond{0, 2, core::BondOrder::SINGLE},
        core::Bond{0, 3, core::BondOrder::SINGLE}, core::Bond{0, 4, core::BondOrder::SINGLE}};

    return core::Molecule{std::move(atoms), std::move(bonds), {}, "methane"};
}

auto make_formally_charged_pair() -> core::Molecule {
    std::vector atoms{core::Atom{7, 1, "N"}, core::Atom{17, -1, "Cl"}};

    return core::Molecule{std::move(atoms), {}, {}, "charged-pair"};
}

auto relabel_atoms(const core::Molecule& molecule, std::span<const std::size_t> new_atom_order)
    -> core::Molecule {
    std::vector<core::Atom> atoms;
    atoms.reserve(molecule.atom_count());

    for (const std::size_t old_index : new_atom_order) {
        atoms.push_back(molecule.atom(old_index));
    }

    std::vector<core::Bond> bonds;
    bonds.reserve(molecule.bond_count());

    std::vector<std::size_t> old_to_new(molecule.atom_count());
    for (std::size_t new_index = 0; new_index < new_atom_order.size(); ++new_index) {
        old_to_new[new_atom_order[new_index]] = new_index;
    }

    for (const auto& bond : molecule.bonds()) {
        bonds.emplace_back(old_to_new[bond.first_atom_index()],
                           old_to_new[bond.second_atom_index()], bond.order());
    }

    std::vector<core::Conformer> conformers;
    conformers.reserve(molecule.conformer_count());

    for (const auto& conformer : molecule.conformers()) {
        std::vector<core::Position> positions;
        positions.reserve(conformer.size());

        for (const std::size_t old_index : new_atom_order) {
            positions.push_back(conformer[old_index]);
        }

        conformers.emplace_back(std::move(positions), std::string{conformer.name()});
    }

    return core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers),
                          std::string{molecule.name()}};
}

auto flip_bond_directions(const core::Molecule& molecule) -> core::Molecule {
    std::vector<core::Bond> bonds;
    bonds.reserve(molecule.bond_count());

    for (const auto& bond : molecule.bonds()) {
        bonds.emplace_back(bond.second_atom_index(), bond.first_atom_index(), bond.order());
    }

    std::vector<core::Atom> atoms{molecule.atoms().begin(), molecule.atoms().end()};
    std::vector<core::Conformer> conformers{molecule.conformers().begin(),
                                            molecule.conformers().end()};

    return core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers),
                          std::string{molecule.name()}};
}

} // namespace chargefw::test
