#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>

#include <print>

using namespace chargefw;

auto make_water() -> core::Molecule {
    std::vector atoms{core::Atom{8, 0, "O"}, core::Atom{1, 0, "H1"}, core::Atom{1, 0, "H2"}};

    std::vector bonds{core::Bond{0, 1, core::BondOrder::SINGLE},
                      core::Bond{0, 2, core::BondOrder::SINGLE}};

    std::vector positions{core::Position{.x = 0.0000, .y = 0.0000, .z = 0.0000},
                                          core::Position{.x = 0.9572, .y = 0.0000, .z = 0.0000},
                                          core::Position{.x = -0.2390, .y = 0.9270, .z = 0.0000}};

    std::vector conformers{core::Conformer{std::move(positions), "conformer-1"}};

    return core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers), "water"};
}

auto main() -> int {

    const auto a = core::Atom{1, 0, "H1"};
    std::println("Atomic number: {}, Formal charge: {}, Name: {}", a.atomic_number(),
                 a.formal_charge(), a.name());

    const auto m = make_water();

    std::println("Name: {}", m.name());

    return 0;
}