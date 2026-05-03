#include "support/test_molecules.h"

#include <chargefw/charges/charge_collection.h>
#include <chargefw/charges/charge_validation.h>
#include <chargefw/core/molecule_collection.h>

#include <cassert>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace charges = chargefw::charges;
namespace core = chargefw::core;

auto main() -> int {
    const charges::ChargeTarget molecule_target{};
    assert(!molecule_target.is_conformer_specific());

    const charges::ChargeTarget conformer_target{.molecule_index = 0, .conformer_index = 0};
    assert(conformer_target.is_conformer_specific());

    const charges::ChargeAssignment water_assignment{
        .target = conformer_target, .charges = charges::AtomicCharges{{-0.8, 0.4, 0.4}}};
    const charges::ChargeAssignment pair_assignment{
        .target = charges::ChargeTarget{.molecule_index = 1, .conformer_index = std::nullopt},
        .charges = charges::AtomicCharges{{1.0, -1.0}}};

    const charges::ChargeSet charge_set{
        "formal", {water_assignment, pair_assignment}, std::string{"default"}};

    assert(charge_set.method_id() == std::string_view{"formal"});
    assert(charge_set.parameter_set_id().has_value());
    assert(*charge_set.parameter_set_id() == std::string_view{"default"});
    assert(charge_set.size() == 2);
    assert(!charge_set.empty());
    assert(charge_set.assignments().size() == 2);
    assert(charge_set.assignment(0).target.molecule_index == 0);
    assert(charge_set.assignment(1).target.molecule_index == 1);

    const charges::ChargeSet parameterless_charge_set{"dummy", {}};
    assert(!parameterless_charge_set.parameter_set_id().has_value());
    assert(parameterless_charge_set.empty());

    const charges::ChargeCollection collection{{charge_set}};
    assert(collection.size() == 1);
    assert(!collection.empty());
    assert(collection.charge_sets().size() == 1);
    assert(collection[0].method_id() == std::string_view{"formal"});

    const charges::ChargeCollection empty_collection{std::vector<charges::ChargeSet>{}};
    assert(empty_collection.empty());

    const core::MoleculeCollection molecules{
        {chargefw::test::make_water(), chargefw::test::make_formally_charged_pair()}, "molecules"};
    charges::validate_charge_collection(molecules, collection);

    bool rejected_empty_method_id = false;

    try {
        [[maybe_unused]] const charges::ChargeSet invalid{"", {}};
    } catch (const std::invalid_argument&) {
        rejected_empty_method_id = true;
    }

    assert(rejected_empty_method_id);

    bool rejected_bad_assignment_index = false;

    try {
        [[maybe_unused]] const auto& invalid = charge_set.assignment(2);
    } catch (const std::out_of_range&) {
        rejected_bad_assignment_index = true;
    }

    assert(rejected_bad_assignment_index);

    bool rejected_bad_charge_set_index = false;

    try {
        [[maybe_unused]] const auto& invalid = collection.at(1);
    } catch (const std::out_of_range&) {
        rejected_bad_charge_set_index = true;
    }

    assert(rejected_bad_charge_set_index);

    bool rejected_bad_molecule_target = false;

    try {
        charges::validate_charge_collection(
            molecules, charges::ChargeCollection{{charges::ChargeSet{
                           "dummy",
                           {{.target = charges::ChargeTarget{.molecule_index = 2},
                             .charges = charges::AtomicCharges{{0.0}}}}}}});
    } catch (const std::invalid_argument&) {
        rejected_bad_molecule_target = true;
    }

    assert(rejected_bad_molecule_target);

    bool rejected_bad_conformer_target = false;

    try {
        charges::validate_charge_collection(
            molecules,
            charges::ChargeCollection{{charges::ChargeSet{
                "dummy",
                {{.target = charges::ChargeTarget{.molecule_index = 1, .conformer_index = 0},
                  .charges = charges::AtomicCharges{{1.0, -1.0}}}}}}});
    } catch (const std::invalid_argument&) {
        rejected_bad_conformer_target = true;
    }

    assert(rejected_bad_conformer_target);

    bool rejected_bad_charge_count = false;

    try {
        charges::validate_charge_collection(
            molecules, charges::ChargeCollection{{charges::ChargeSet{
                           "dummy",
                           {{.target = charges::ChargeTarget{.molecule_index = 0},
                             .charges = charges::AtomicCharges{{0.0, 0.0}}}}}}});
    } catch (const std::invalid_argument&) {
        rejected_bad_charge_count = true;
    }

    assert(rejected_bad_charge_count);

    return 0;
}
