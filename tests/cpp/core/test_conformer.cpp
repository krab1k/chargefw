#include <chargefw/core/conformer.h>
#include <chargefw/core/position.h>

#include <cassert>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace core = chargefw::core;

auto main() -> int
{
    const core::Conformer conformer{
        {
            core::Position{0.0, 0.0, 0.0},
            core::Position{1.0, 0.0, 0.0}
        },
        "model-1"
    };

    assert(conformer.name() == std::string_view{"model-1"});
    assert(conformer.size() == 2);
    assert(!conformer.empty());
    assert(conformer.positions().size() == 2);
    assert(conformer[1].x == 1.0);
    assert(conformer.at(0).z == 0.0);

    const core::Conformer empty_conformer{std::vector<core::Position>{}, "empty"};
    assert(empty_conformer.empty());
    assert(empty_conformer.size() == 0);
    assert(empty_conformer.positions().empty());

    bool rejected_bad_index = false;

    try {
        [[maybe_unused]] const auto& invalid = conformer.at(2);
    } catch (const std::out_of_range&) {
        rejected_bad_index = true;
    }

    assert(rejected_bad_index);

    return 0;
}
