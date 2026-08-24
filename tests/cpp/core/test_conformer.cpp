#include <chargefw/core/conformer.h>
#include <chargefw/core/position.h>

#include <stdexcept>
#include <string_view>
#include <vector>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;

TEST_CASE("conformer exposes named positions", "[core][conformer]") {
    const core::Conformer conformer{{core::Position{.x = 0.0, .y = 0.0, .z = 0.0},
                                     core::Position{.x = 1.0, .y = 0.0, .z = 0.0}},
                                    "model-1"};

    CHECK(conformer.name() == std::string_view{"model-1"});
    CHECK(conformer.size() == 2);
    CHECK_FALSE(conformer.empty());
    CHECK(conformer.positions().size() == 2);
    CHECK(conformer[1].x == 1.0);
    CHECK(conformer.at(0).z == 0.0);
}

TEST_CASE("conformer supports an empty position list", "[core][conformer]") {
    const core::Conformer conformer{std::vector<core::Position>{}, "empty"};

    CHECK(conformer.empty());
    CHECK(conformer.positions().empty());
}

TEST_CASE("conformer bounds-checked access rejects an invalid index", "[core][conformer]") {
    const core::Conformer conformer{{core::Position{.x = 0.0, .y = 0.0, .z = 0.0}}, "model-1"};

    CHECK_THROWS_AS(conformer.at(1), std::out_of_range);
}
