#include "support/test_assertions.h"

#include <chargefw/calculation/execution_policy.h>

#include <cassert>
#include <limits>
#include <optional>

namespace calculation = chargefw::calculation;

auto main() -> int {
    assert(calculation::execution_selection_kind_from_string("auto") ==
           calculation::ExecutionSelectionKind::automatic);
    assert(calculation::execution_selection_kind_from_string("full") ==
           calculation::ExecutionSelectionKind::full);
    assert(calculation::execution_selection_kind_from_string("cutoff") ==
           calculation::ExecutionSelectionKind::cutoff);
    assert(calculation::execution_selection_kind_from_string("cover") ==
           calculation::ExecutionSelectionKind::cover);
    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(calculation::execution_selection_kind_from_string("unknown"));
    }));
    assert(calculation::charge_correction_policy_from_string("none") ==
           calculation::ChargeCorrectionPolicy::none);
    assert(calculation::charge_correction_policy_from_string("uniform") ==
           calculation::ChargeCorrectionPolicy::uniform);
    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(calculation::charge_correction_policy_from_string("unknown"));
    }));
    assert(calculation::to_string(calculation::ExecutionSelectionKind::automatic) == "auto");
    assert(calculation::to_string(calculation::ExecutionSelectionKind::full) == "full");
    assert(calculation::to_string(calculation::ExecutionMode::cutoff) == "cutoff");
    assert(calculation::to_string(calculation::ChargeCorrectionPolicy::uniform) == "uniform");
    const calculation::ExecutionPolicy default_policy;
    assert(default_policy.mode() == calculation::ExecutionMode::full);
    assert(!default_policy.radius().has_value());

    const calculation::ExecutionPolicy full_policy{calculation::ExecutionMode::full};
    assert(full_policy.mode() == calculation::ExecutionMode::full);
    assert(!full_policy.radius().has_value());
    assert(full_policy.charge_correction() == calculation::ChargeCorrectionPolicy::none);
    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::full, 8.0});
    }));

    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::cutoff});
    }));
    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::cover,
                                                       std::numeric_limits<double>::quiet_NaN()});
    }));
    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::cutoff,
                                                       std::numeric_limits<double>::infinity()});
    }));
    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::cover, 7.99});
    }));

    const calculation::ExecutionPolicy cutoff_policy{calculation::ExecutionMode::cutoff, 8.0};
    assert(cutoff_policy.radius() == std::optional<double>{8.0});
    assert(cutoff_policy.charge_correction() == calculation::ChargeCorrectionPolicy::none);
    const calculation::ExecutionPolicy corrected_cutoff_policy{
        calculation::ExecutionMode::cutoff, 8.0, calculation::ChargeCorrectionPolicy::uniform};
    assert(corrected_cutoff_policy.charge_correction() ==
           calculation::ChargeCorrectionPolicy::uniform);
    const calculation::ExecutionPolicy cover_policy{calculation::ExecutionMode::cover, 12.0};
    assert(cover_policy.radius() == std::optional<double>{12.0});

    const calculation::ExecutionSelection default_selection;
    assert(default_selection.kind() == calculation::ExecutionSelectionKind::automatic);
    assert(!default_selection.radius().has_value());
    const calculation::ExecutionSelection automatic_radius{
        calculation::ExecutionSelectionKind::automatic, 8.0};
    assert(automatic_radius.radius() == std::optional<double>{8.0});
    const calculation::ExecutionSelection full_selection{calculation::ExecutionSelectionKind::full};
    assert(!full_selection.radius().has_value());
    assert(!full_selection.charge_correction().has_value());
    const calculation::ExecutionSelection corrected_cutoff_selection{
        calculation::ExecutionSelectionKind::cutoff, 8.0,
        calculation::ChargeCorrectionPolicy::uniform};
    assert(corrected_cutoff_selection.charge_correction() ==
           std::optional{calculation::ChargeCorrectionPolicy::uniform});
    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::automatic, 7.99});
    }));
    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full, 8.0});
    }));
    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff});
    }));
    assert(chargefw::test::throws_invalid_argument([] -> void {
        static_cast<void>(
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full, std::nullopt,
                                            calculation::ChargeCorrectionPolicy::uniform});
    }));

    const calculation::ResourcePolicy default_resources;
    assert(default_resources.cutoff_atom_threshold ==
           std::optional<std::size_t>{calculation::default_cutoff_atom_threshold});
    assert(default_resources.cover_atom_threshold ==
           std::optional<std::size_t>{calculation::default_cover_atom_threshold});
    const calculation::ResourcePolicy finite_resources{.cutoff_atom_threshold = 42,
                                                       .cover_atom_threshold = 84};
    assert(finite_resources.cutoff_atom_threshold == std::optional<std::size_t>{42});
    assert(finite_resources.cover_atom_threshold == std::optional<std::size_t>{84});
    const calculation::ResourcePolicy unlimited_resources{.cutoff_atom_threshold = std::nullopt,
                                                          .cover_atom_threshold = std::nullopt};
    assert(!unlimited_resources.cutoff_atom_threshold.has_value());
    assert(!unlimited_resources.cover_atom_threshold.has_value());
    assert(calculation::default_automatic_reduced_radius == 12.0);

    return 0;
}
