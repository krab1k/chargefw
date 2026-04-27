#include <chargefw/core/atom.h>

#include <cassert>
#include <stdexcept>
#include <string_view>

auto main() -> int {
    const chargefw::core::Atom carbon{6, 0, "C1"};

    assert(carbon.atomic_number() == 6);
    assert(carbon.formal_charge() == 0);
    assert(carbon.name() == std::string_view{"C1"});

    bool threw = false;

    try {
        [[maybe_unused]] const chargefw::core::Atom invalid{0};
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    assert(threw);

    return 0;
}