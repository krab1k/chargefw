#include <chargefw/calculation/execution_policy.h>

#include <limits>
#include <optional>
#include <stdexcept>

#include <snitch/snitch.hpp>

namespace calculation = chargefw::calculation;

TEST_CASE("execution selection and charge correction parse from strings",
          "[calculation][execution-policy]") {
    CHECK(calculation::execution_selection_kind_from_string("auto") ==
          calculation::ExecutionSelectionKind::automatic);
    CHECK(calculation::execution_selection_kind_from_string("full") ==
          calculation::ExecutionSelectionKind::full);
    CHECK(calculation::execution_selection_kind_from_string("cutoff") ==
          calculation::ExecutionSelectionKind::cutoff);
    CHECK(calculation::execution_selection_kind_from_string("cover") ==
          calculation::ExecutionSelectionKind::cover);

    const auto bad_selection = [] {
        static_cast<void>(calculation::execution_selection_kind_from_string("unknown"));
    };
    CHECK_THROWS_AS(bad_selection(), std::invalid_argument);

    CHECK(calculation::charge_correction_policy_from_string("none") ==
          calculation::ChargeCorrectionPolicy::none);
    CHECK(calculation::charge_correction_policy_from_string("uniform") ==
          calculation::ChargeCorrectionPolicy::uniform);

    const auto bad_correction = [] {
        static_cast<void>(calculation::charge_correction_policy_from_string("unknown"));
    };
    CHECK_THROWS_AS(bad_correction(), std::invalid_argument);
}

TEST_CASE("execution policy and selection convert to strings", "[calculation][execution-policy]") {
    CHECK(calculation::to_string(calculation::ExecutionSelectionKind::automatic) == "auto");
    CHECK(calculation::to_string(calculation::ExecutionSelectionKind::full) == "full");
    CHECK(calculation::to_string(calculation::ExecutionMode::cutoff) == "cutoff");
    CHECK(calculation::to_string(calculation::ChargeCorrectionPolicy::uniform) == "uniform");
}

TEST_CASE("execution policy validates mode, radius, and correction",
          "[calculation][execution-policy]") {
    const calculation::ExecutionPolicy default_policy;
    CHECK(default_policy.mode() == calculation::ExecutionMode::full);
    CHECK_FALSE(default_policy.radius().has_value());

    const calculation::ExecutionPolicy full_policy{calculation::ExecutionMode::full};
    CHECK(full_policy.mode() == calculation::ExecutionMode::full);
    CHECK_FALSE(full_policy.radius().has_value());
    CHECK(full_policy.charge_correction() == calculation::ChargeCorrectionPolicy::none);

    const auto full_with_radius = [] {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::full, 8.0});
    };
    CHECK_THROWS_AS(full_with_radius(), std::invalid_argument);

    const auto cutoff_no_radius = [] {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::cutoff});
    };
    CHECK_THROWS_AS(cutoff_no_radius(), std::invalid_argument);

    const auto cutoff_nan = [] {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::cutoff,
                                                       std::numeric_limits<double>::quiet_NaN()});
    };
    CHECK_THROWS_AS(cutoff_nan(), std::invalid_argument);

    const auto cutoff_infinite = [] {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::cutoff,
                                                       std::numeric_limits<double>::infinity()});
    };
    CHECK_THROWS_AS(cutoff_infinite(), std::invalid_argument);

    const auto cover_too_small = [] {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::cover, 7.99});
    };
    CHECK_THROWS_AS(cover_too_small(), std::invalid_argument);

    const calculation::ExecutionPolicy cutoff_policy{calculation::ExecutionMode::cutoff, 8.0};
    CHECK(cutoff_policy.radius() == std::optional<double>{8.0});
    CHECK(cutoff_policy.charge_correction() == calculation::ChargeCorrectionPolicy::none);

    const calculation::ExecutionPolicy corrected_cutoff_policy{
        calculation::ExecutionMode::cutoff, 8.0, calculation::ChargeCorrectionPolicy::uniform};
    CHECK(corrected_cutoff_policy.charge_correction() ==
          calculation::ChargeCorrectionPolicy::uniform);

    const calculation::ExecutionPolicy cover_policy{calculation::ExecutionMode::cover, 12.0};
    CHECK(cover_policy.radius() == std::optional<double>{12.0});
}

TEST_CASE("execution selection validates kind, radius, and correction",
          "[calculation][execution-policy]") {
    const calculation::ExecutionSelection default_selection;
    CHECK(default_selection.kind() == calculation::ExecutionSelectionKind::automatic);
    CHECK_FALSE(default_selection.radius().has_value());

    const calculation::ExecutionSelection automatic_radius{
        calculation::ExecutionSelectionKind::automatic, 8.0};
    CHECK(automatic_radius.radius() == std::optional<double>{8.0});

    const calculation::ExecutionSelection full_selection{calculation::ExecutionSelectionKind::full};
    CHECK_FALSE(full_selection.radius().has_value());
    CHECK_FALSE(full_selection.charge_correction().has_value());

    const calculation::ExecutionSelection corrected_cutoff_selection{
        calculation::ExecutionSelectionKind::cutoff, 8.0,
        calculation::ChargeCorrectionPolicy::uniform};
    CHECK(corrected_cutoff_selection.charge_correction() ==
          std::optional{calculation::ChargeCorrectionPolicy::uniform});

    const auto auto_too_small = [] {
        static_cast<void>(
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::automatic, 7.99});
    };
    CHECK_THROWS_AS(auto_too_small(), std::invalid_argument);

    const auto full_with_radius = [] {
        static_cast<void>(
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full, 8.0});
    };
    CHECK_THROWS_AS(full_with_radius(), std::invalid_argument);

    const auto cutoff_no_radius = [] {
        static_cast<void>(
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff});
    };
    CHECK_THROWS_AS(cutoff_no_radius(), std::invalid_argument);

    const auto full_with_correction = [] {
        static_cast<void>(
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full, std::nullopt,
                                            calculation::ChargeCorrectionPolicy::uniform});
    };
    CHECK_THROWS_AS(full_with_correction(), std::invalid_argument);
}

TEST_CASE("resource policy exposes thresholds", "[calculation][execution-policy]") {
    const calculation::ResourcePolicy default_resources;
    CHECK(default_resources.cutoff_atom_threshold ==
          std::optional<std::size_t>{calculation::default_cutoff_atom_threshold});
    CHECK(default_resources.cover_atom_threshold ==
          std::optional<std::size_t>{calculation::default_cover_atom_threshold});

    const calculation::ResourcePolicy finite_resources{.cutoff_atom_threshold = 42,
                                                       .cover_atom_threshold = 84};
    CHECK(finite_resources.cutoff_atom_threshold == std::optional<std::size_t>{42});
    CHECK(finite_resources.cover_atom_threshold == std::optional<std::size_t>{84});

    const calculation::ResourcePolicy unlimited_resources{.cutoff_atom_threshold = std::nullopt,
                                                          .cover_atom_threshold = std::nullopt};
    CHECK_FALSE(unlimited_resources.cutoff_atom_threshold.has_value());
    CHECK_FALSE(unlimited_resources.cover_atom_threshold.has_value());

    CHECK(calculation::default_automatic_reduced_radius == 12.0);
}
