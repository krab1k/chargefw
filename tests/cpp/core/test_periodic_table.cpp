#include <chargefw/core/periodic_table.h>

#include <optional>
#include <stdexcept>
#include <string_view>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;

TEST_CASE("periodic table looks up known elements", "[core][periodic-table]") {
    const auto& table = core::periodic_table();

    CHECK(table.contains(1));
    CHECK(table.contains("H"));
    CHECK_FALSE(table.contains(0));
    CHECK(table.contains(100));
    CHECK_FALSE(table.contains(101));
    CHECK_FALSE(table.contains("Xx"));
    CHECK(table.element(1).symbol == std::string_view{"H"});
    CHECK(table.element("C").atomic_number == 6);
    CHECK(core::element_symbol(8) == std::string_view{"O"});
}

TEST_CASE("periodic table reports main-group valence electrons", "[core][periodic-table]") {
    const auto& table = core::periodic_table();

    CHECK(table.element("H").valence_electron_count() == std::optional{1});
    CHECK(table.element("He").valence_electron_count() == std::optional{2});
    CHECK(table.element("C").valence_electron_count() == std::optional{4});
    CHECK(table.element("N").valence_electron_count() == std::optional{5});
    CHECK(table.element("O").valence_electron_count() == std::optional{6});
    CHECK(table.element("F").valence_electron_count() == std::optional{7});
    CHECK(table.element("Ne").valence_electron_count() == std::optional{8});
    CHECK_FALSE(table.element("Fe").valence_electron_count().has_value());
    CHECK_FALSE(table.element("Zn").valence_electron_count().has_value());
}

TEST_CASE("periodic table rejects unknown elements", "[core][periodic-table]") {
    const auto& table = core::periodic_table();

    CHECK_THROWS_AS(table.element(0), std::out_of_range);
    CHECK(table.element(100).symbol == std::string_view{"Fm"});
    CHECK_THROWS_AS(table.element(101), std::out_of_range);
    CHECK_THROWS_AS(table.element("Xx"), std::out_of_range);
}
