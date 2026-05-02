#include <chargefw/core/atom.h>

#include <cassert>
#include <stdexcept>
#include <string_view>

namespace core = chargefw::core;

auto main() -> int {
    const core::Atom carbon{6, 0, "C1"};

    assert(carbon.atomic_number() == 6);
    assert(carbon.formal_charge() == 0);
    assert(carbon.name() == std::string_view{"C1"});

    const core::Atom ammonium_nitrogen{7, 1, "N"};

    assert(ammonium_nitrogen.atomic_number() == 7);
    assert(ammonium_nitrogen.formal_charge() == 1);
    assert(ammonium_nitrogen.name() == std::string_view{"N"});

    const core::Atom unnamed_oxygen{8};

    assert(unnamed_oxygen.atomic_number() == 8);
    assert(unnamed_oxygen.formal_charge() == 0);
    assert(unnamed_oxygen.name().empty());

    bool rejected_zero = false;

    try {
        [[maybe_unused]] const core::Atom invalid{0};
    } catch (const std::invalid_argument&) {
        rejected_zero = true;
    }

    assert(rejected_zero);

    bool rejected_too_large = false;

    try {
        [[maybe_unused]] const core::Atom invalid{119};
    } catch (const std::invalid_argument&) {
        rejected_too_large = true;
    }

    assert(rejected_too_large);

    return 0;
}