#include <chargefw/core/molecule.h>

#include <vector>

auto main() -> int {
    const chargefw::core::Molecule molecule{std::vector<chargefw::core::Atom>{}};
    return molecule.atom_count() == 0 ? 0 : 1;
}
