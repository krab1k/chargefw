#include <chargefw/calculation/execution_policy.h>

#include <cassert>
#include <limits>
#include <optional>
#include <stdexcept>

namespace calculation = chargefw::calculation;

namespace {

template <typename Callable> auto throws_invalid_argument(Callable&& callable) -> bool {
    try {
        callable();
    } catch (const std::invalid_argument&) {
        return true;
    }

    return false;
}

} // namespace

auto main() -> int {
    const calculation::ExecutionPolicy default_policy;
    assert(default_policy.mode() == calculation::ExecutionMode::full);
    assert(!default_policy.radius().has_value());

    const calculation::ExecutionPolicy full_policy{calculation::ExecutionMode::full};
    assert(full_policy.mode() == calculation::ExecutionMode::full);
    assert(!full_policy.radius().has_value());
    assert(throws_invalid_argument([] -> void {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::full, 8.0});
    }));

    assert(throws_invalid_argument([] -> void {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::cutoff});
    }));
    assert(throws_invalid_argument([] -> void {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::cover,
                                                       std::numeric_limits<double>::quiet_NaN()});
    }));
    assert(throws_invalid_argument([] -> void {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::cutoff,
                                                       std::numeric_limits<double>::infinity()});
    }));
    assert(throws_invalid_argument([] -> void {
        static_cast<void>(calculation::ExecutionPolicy{calculation::ExecutionMode::cover, 7.99});
    }));

    const calculation::ExecutionPolicy cutoff_policy{calculation::ExecutionMode::cutoff, 8.0};
    assert(cutoff_policy.radius() == std::optional<double>{8.0});
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
    assert(throws_invalid_argument([] -> void {
        static_cast<void>(
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::automatic, 7.99});
    }));
    assert(throws_invalid_argument([] -> void {
        static_cast<void>(
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full, 8.0});
    }));
    assert(throws_invalid_argument([] -> void {
        static_cast<void>(
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff});
    }));

    const calculation::ResourcePolicy default_resources;
    assert(default_resources.full_atom_threshold ==
           std::optional<std::size_t>{calculation::default_full_atom_threshold});
    const calculation::ResourcePolicy finite_resources{.full_atom_threshold = 42};
    assert(finite_resources.full_atom_threshold == std::optional<std::size_t>{42});
    const calculation::ResourcePolicy unlimited_resources{.full_atom_threshold = std::nullopt};
    assert(!unlimited_resources.full_atom_threshold.has_value());

    return 0;
}
