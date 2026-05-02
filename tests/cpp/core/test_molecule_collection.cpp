#include "support/test_molecules.h"

#include <chargefw/core/molecule_collection.h>

#include <cassert>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace core = chargefw::core;

auto main() -> int {
    const core::MoleculeCollection collection{
        {chargefw::test::make_water(), chargefw::test::make_formally_charged_pair()}, "examples"};

    assert(collection.name() == std::string_view{"examples"});
    assert(collection.size() == 2);
    assert(!collection.empty());
    assert(collection.molecules().size() == 2);
    assert(collection[0].name() == std::string_view{"water"});
    assert(collection.at(1).name() == std::string_view{"charged-pair"});

    const core::MoleculeCollection empty_collection{std::vector<core::Molecule>{}, "empty"};
    assert(empty_collection.empty());
    assert(empty_collection.size() == 0);

    bool rejected_bad_index = false;

    try {
        [[maybe_unused]] const auto& invalid = collection.at(2);
    } catch (const std::out_of_range&) {
        rejected_bad_index = true;
    }

    assert(rejected_bad_index);

    return 0;
}
