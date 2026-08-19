#include "support/test_molecules.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;

namespace {

auto approximately_equal(const double first, const double second, const double tolerance) -> bool {
    return std::abs(first - second) <= tolerance;
}

} // namespace

auto main() -> int {
    static_assert(!std::is_constructible_v<features::ConformerFeatures, core::Molecule&&>);
    static_assert(
        !std::is_constructible_v<features::ConformerFeatures, core::Molecule&&, std::size_t>);

    const auto water = chargefw::test::make_water();
    const features::ConformerFeatures geometry{water, 0};

    assert(&geometry.molecule() == &water);
    assert(geometry.conformer_index() == 0);

    const auto& oxygen = geometry.position(0);
    assert(approximately_equal(oxygen.x, 0.0, 1e-12));
    assert(approximately_equal(oxygen.y, 0.0, 1e-12));
    assert(approximately_equal(oxygen.z, 0.0, 1e-12));

    assert(approximately_equal(geometry.squared_distance(0, 1), 0.9572 * 0.9572, 1e-12));
    assert(approximately_equal(geometry.distance(0, 1), 0.9572, 1e-12));

    bool rejected_bad_conformer = false;

    try {
        [[maybe_unused]] const features::ConformerFeatures invalid{water, 1};
    } catch (const std::out_of_range&) {
        rejected_bad_conformer = true;
    }

    assert(rejected_bad_conformer);

    bool rejected_bad_atom = false;

    try {
        [[maybe_unused]] const auto& invalid_position = geometry.position(3);
    } catch (const std::out_of_range&) {
        rejected_bad_atom = true;
    }

    assert(rejected_bad_atom);

    return 0;
}
