#include <chargefw/core/periodic_table.h>

#include <cassert>
#include <optional>
#include <stdexcept>
#include <string_view>

namespace core = chargefw::core;

auto main() -> int {
    const auto& table = core::periodic_table();

    assert(table.contains(1));
    assert(table.contains("H"));
    assert(!table.contains(0));
    assert(!table.contains("Xx"));

    assert(table.element(1).symbol == std::string_view{"H"});
    assert(table.element("C").atomic_number == 6);
    assert(core::element_symbol(8) == std::string_view{"O"});

    assert(table.element("H").valence_electron_count() == std::optional{1});
    assert(table.element("He").valence_electron_count() == std::optional{2});
    assert(table.element("C").valence_electron_count() == std::optional{4});
    assert(table.element("N").valence_electron_count() == std::optional{5});
    assert(table.element("O").valence_electron_count() == std::optional{6});
    assert(table.element("F").valence_electron_count() == std::optional{7});
    assert(table.element("Ne").valence_electron_count() == std::optional{8});

    assert(!table.element("Fe").valence_electron_count().has_value());
    assert(!table.element("Zn").valence_electron_count().has_value());

    bool rejected_unknown_atomic_number = false;

    try {
        [[maybe_unused]] const auto& invalid = table.element(0);
    } catch (const std::out_of_range&) {
        rejected_unknown_atomic_number = true;
    }

    assert(rejected_unknown_atomic_number);

    bool rejected_unknown_symbol = false;

    try {
        [[maybe_unused]] const auto& invalid = table.element("Xx");
    } catch (const std::out_of_range&) {
        rejected_unknown_symbol = true;
    }

    assert(rejected_unknown_symbol);

    return 0;
}
