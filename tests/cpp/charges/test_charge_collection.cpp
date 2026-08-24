#include "support/test_molecules.h"

#include <chargefw/charges/charge_collection.h>
#include <chargefw/charges/charge_validation.h>
#include <chargefw/core/molecule_collection.h>

#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <snitch/snitch.hpp>

namespace charges = chargefw::charges;
namespace core = chargefw::core;

TEST_CASE("charge target distinguishes conformer-specific and molecule-level",
          "[charges][charge-collection]") {
    const charges::ChargeTarget molecule_target{};

    CHECK_FALSE(molecule_target.is_conformer_specific());

    const charges::ChargeTarget conformer_target{.molecule_index = 0, .conformer_index = 0};

    CHECK(conformer_target.is_conformer_specific());
}

TEST_CASE("charge set stores method, parameter, and assignments", "[charges][charge-collection]") {
    const charges::ChargeTarget conformer_target{.molecule_index = 0, .conformer_index = 0};
    const charges::ChargeAssignment water_assignment{
        .target = conformer_target, .charges = charges::AtomicCharges{{-0.8, 0.4, 0.4}}};
    const charges::ChargeAssignment pair_assignment{
        .target = charges::ChargeTarget{.molecule_index = 1, .conformer_index = std::nullopt},
        .charges = charges::AtomicCharges{{1.0, -1.0}}};

    const charges::ChargeSet charge_set{
        "formal", {water_assignment, pair_assignment}, std::string{"default"}};

    CHECK(charge_set.method_id() == std::string_view{"formal"});
    REQUIRE(charge_set.parameter_set_id().has_value());
    CHECK(*charge_set.parameter_set_id() == std::string_view{"default"});
    REQUIRE(charge_set.size() == 2);
    CHECK_FALSE(charge_set.empty());
    REQUIRE(charge_set.assignments().size() == 2);
    CHECK(charge_set.assignment(0).target.molecule_index == 0);
    CHECK(charge_set.assignment(1).target.molecule_index == 1);
}

TEST_CASE("parameterless charge set has no parameter-set id", "[charges][charge-collection]") {
    const charges::ChargeSet parameterless_charge_set{"dummy", {}};

    CHECK_FALSE(parameterless_charge_set.parameter_set_id().has_value());
    CHECK(parameterless_charge_set.empty());
}

TEST_CASE("charge collection preserves ordered charge sets", "[charges][charge-collection]") {
    const charges::ChargeTarget conformer_target{.molecule_index = 0, .conformer_index = 0};
    const charges::ChargeAssignment water_assignment{
        .target = conformer_target, .charges = charges::AtomicCharges{{-0.8, 0.4, 0.4}}};
    const charges::ChargeAssignment pair_assignment{
        .target = charges::ChargeTarget{.molecule_index = 1, .conformer_index = std::nullopt},
        .charges = charges::AtomicCharges{{1.0, -1.0}}};
    const charges::ChargeSet charge_set{
        "formal", {water_assignment, pair_assignment}, std::string{"default"}};
    const charges::ChargeCollection collection{{charge_set}};

    REQUIRE(collection.size() == 1);
    CHECK_FALSE(collection.empty());
    REQUIRE(collection.charge_sets().size() == 1);
    CHECK(collection[0].method_id() == std::string_view{"formal"});

    const charges::ChargeCollection empty_collection{std::vector<charges::ChargeSet>{}};
    CHECK(empty_collection.empty());
}

TEST_CASE("charge collection validates against molecule collection",
          "[charges][charge-collection]") {
    const charges::ChargeTarget conformer_target{.molecule_index = 0, .conformer_index = 0};
    const charges::ChargeAssignment water_assignment{
        .target = conformer_target, .charges = charges::AtomicCharges{{-0.8, 0.4, 0.4}}};
    const charges::ChargeAssignment pair_assignment{
        .target = charges::ChargeTarget{.molecule_index = 1, .conformer_index = std::nullopt},
        .charges = charges::AtomicCharges{{1.0, -1.0}}};
    const charges::ChargeSet charge_set{
        "formal", {water_assignment, pair_assignment}, std::string{"default"}};
    const charges::ChargeCollection collection{{charge_set}};
    const core::MoleculeCollection molecules{
        {chargefw::test::make_water(), chargefw::test::make_formally_charged_pair()}, "molecules"};

    CHECK_NOTHROW(charges::validate_charge_collection(molecules, collection));
}

TEST_CASE("charge set and collection reject empty method id and invalid indices",
          "[charges][charge-collection]") {
    const charges::ChargeTarget conformer_target{.molecule_index = 0, .conformer_index = 0};
    const charges::ChargeAssignment water_assignment{
        .target = conformer_target, .charges = charges::AtomicCharges{{-0.8, 0.4, 0.4}}};
    const charges::ChargeAssignment pair_assignment{
        .target = charges::ChargeTarget{.molecule_index = 1, .conformer_index = std::nullopt},
        .charges = charges::AtomicCharges{{1.0, -1.0}}};
    const charges::ChargeSet charge_set{
        "formal", {water_assignment, pair_assignment}, std::string{"default"}};
    const charges::ChargeCollection collection{{charge_set}};

    CHECK_THROWS_AS((charges::ChargeSet{"", {}}), std::invalid_argument);
    CHECK_THROWS_AS(charge_set.assignment(2), std::out_of_range);
    CHECK_THROWS_AS(collection.at(1), std::out_of_range);
}

TEST_CASE("charge validation rejects bad targets and charge counts",
          "[charges][charge-collection]") {
    const core::MoleculeCollection molecules{
        {chargefw::test::make_water(), chargefw::test::make_formally_charged_pair()}, "molecules"};

    const auto make_collection = [](std::size_t molecule_index,
                                    std::optional<std::size_t> conformer_index,
                                    std::vector<double> charge_values) {
        return charges::ChargeCollection{{charges::ChargeSet{
            "dummy",
            {{.target = charges::ChargeTarget{.molecule_index = molecule_index,
                                              .conformer_index = conformer_index},
              .charges = charges::AtomicCharges{std::move(charge_values)}}}}}};
    };

    // molecule index out of range
    CHECK_THROWS_AS(
        charges::validate_charge_collection(molecules, make_collection(2, std::nullopt, {0.0})),
        std::invalid_argument);

    // conformer index on a molecule without conformers
    CHECK_THROWS_AS(
        charges::validate_charge_collection(molecules, make_collection(1, 0, {1.0, -1.0})),
        std::invalid_argument);

    // charge count does not match atom count
    CHECK_THROWS_AS(
        charges::validate_charge_collection(molecules, make_collection(0, 0, {0.0, 0.0})),
        std::invalid_argument);
}

TEST_CASE("charge validation rejects out-of-order targets", "[charges][charge-collection]") {
    const core::MoleculeCollection molecules{
        {chargefw::test::make_water(), chargefw::test::make_formally_charged_pair()}, "molecules"};

    const auto make_out_of_order = [] {
        return charges::ChargeCollection{
            {charges::ChargeSet{"dummy",
                                {{.target = charges::ChargeTarget{.molecule_index = 1},
                                  .charges = charges::AtomicCharges{{1.0, -1.0}}},
                                 {.target = charges::ChargeTarget{.molecule_index = 0},
                                  .charges = charges::AtomicCharges{{-0.8, 0.4, 0.4}}}}}}};
    };

    CHECK_THROWS_AS(charges::validate_charge_collection(molecules, make_out_of_order()),
                    std::invalid_argument);
}
