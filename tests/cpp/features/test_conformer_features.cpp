#include "support/test_molecules.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;
namespace features = chargefw::features;

static_assert(!std::is_constructible_v<features::ConformerFeatures, core::Molecule&&>);
static_assert(!std::is_constructible_v<features::ConformerFeatures, core::Molecule&&, std::size_t>);

TEST_CASE("conformer features expose molecule reference and index", "[features][conformer]") {
    const auto water = chargefw::test::make_water();
    const features::ConformerFeatures geometry{water, 0};

    CHECK(&geometry.molecule() == &water);
    CHECK(geometry.conformer_index() == 0);
}

TEST_CASE("conformer features expose positions and distances", "[features][conformer]") {
    const auto water = chargefw::test::make_water();
    const features::ConformerFeatures geometry{water, 0};

    const auto& oxygen = geometry.position(0);
    CHECK(std::abs(oxygen.x - 0.0) < 1e-12);
    CHECK(std::abs(oxygen.y - 0.0) < 1e-12);
    CHECK(std::abs(oxygen.z - 0.0) < 1e-12);

    constexpr auto squared_oh_distance = 0.7570 * 0.7570 + 0.5859 * 0.5859;
    CHECK(std::abs(geometry.squared_distance(0, 1) - squared_oh_distance) < 1e-12);
    CHECK(std::abs(geometry.distance(0, 1) - std::sqrt(squared_oh_distance)) < 1e-12);
    CHECK(std::abs(geometry.distance(0, 1) - geometry.distance(0, 2)) < 1e-12);
}

TEST_CASE("conformer features reject invalid conformer and atom indices", "[features][conformer]") {
    const auto water = chargefw::test::make_water();

    CHECK_THROWS_AS((features::ConformerFeatures{water, 1}), std::out_of_range);

    const features::ConformerFeatures geometry{water, 0};
    CHECK_THROWS_AS(geometry.position(3), std::out_of_range);
}

TEST_CASE("conformer features detect coincident and non-finite coordinates",
          "[features][conformer]") {
    const core::Molecule geometric_edge_cases{
        std::vector{core::Atom{1}, core::Atom{1}},
        {},
        std::vector{core::Conformer{{core::Position{}, core::Position{-0.0}}, "duplicate"},
                    core::Conformer{{core::Position{},
                                     core::Position{.x = std::numeric_limits<double>::infinity()}},
                                    "nonfinite"}}};

    const features::ConformerFeatures duplicate{geometric_edge_cases, 0};
    CHECK_FALSE(duplicate.first_nonfinite_atom_index().has_value());
    const auto coincident = duplicate.coincident_atom_indices();
    REQUIRE(coincident.has_value());
    CHECK(coincident->first == 0);
    CHECK(coincident->second == 1);

    const features::ConformerFeatures nonfinite{geometric_edge_cases, 1};
    CHECK(nonfinite.first_nonfinite_atom_index() == 1);
    CHECK_FALSE(nonfinite.coincident_atom_indices().has_value());
}
