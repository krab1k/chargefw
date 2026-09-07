#include <chargefw/calculation/execution_policy.h>

#include <cmath>
#include <stdexcept>
#include <string>

namespace chargefw::calculation {
namespace {

auto validate_reduced_radius(const std::optional<double> radius) -> void {
    if (!radius.has_value()) {
        throw std::invalid_argument{"cutoff and cover execution require a radius"};
    }

    if (!std::isfinite(*radius)) {
        throw std::invalid_argument{"cutoff and cover radius must be finite"};
    }

    if (*radius < minimum_reduced_radius) {
        throw std::invalid_argument{"cutoff and cover radius must be at least " +
                                    std::to_string(minimum_reduced_radius) + " angstrom"};
    }
}

auto validate_policy(const ExecutionMode mode, const std::optional<double> radius) -> void {
    switch (mode) {
    case ExecutionMode::full:
        if (radius.has_value()) {
            throw std::invalid_argument{"full execution does not accept a radius"};
        }
        return;
    case ExecutionMode::cutoff:
    case ExecutionMode::cover:
        validate_reduced_radius(radius);
        return;
    }

    throw std::invalid_argument{"unknown execution mode"};
}

auto validate_selection(const ExecutionSelectionKind kind, const std::optional<double> radius)
    -> void {
    switch (kind) {
    case ExecutionSelectionKind::automatic:
        if (radius.has_value()) {
            validate_reduced_radius(radius);
        }
        return;
    case ExecutionSelectionKind::full:
        validate_policy(ExecutionMode::full, radius);
        return;
    case ExecutionSelectionKind::cutoff:
    case ExecutionSelectionKind::cover:
        validate_reduced_radius(radius);
        return;
    }

    throw std::invalid_argument{"unknown execution selection"};
}

} // namespace

auto execution_selection_kind_from_string(const std::string_view value) -> ExecutionSelectionKind {
    if (value == "auto") {
        return ExecutionSelectionKind::automatic;
    }
    if (value == "full") {
        return ExecutionSelectionKind::full;
    }
    if (value == "cutoff") {
        return ExecutionSelectionKind::cutoff;
    }
    if (value == "cover") {
        return ExecutionSelectionKind::cover;
    }
    throw std::invalid_argument{"unknown execution selection: " + std::string{value}};
}

auto to_string(const ExecutionSelectionKind value) -> std::string_view {
    switch (value) {
    case ExecutionSelectionKind::automatic:
        return "auto";
    case ExecutionSelectionKind::full:
        return "full";
    case ExecutionSelectionKind::cutoff:
        return "cutoff";
    case ExecutionSelectionKind::cover:
        return "cover";
    }
    throw std::invalid_argument{"unknown execution selection"};
}

auto to_string(const ExecutionMode value) -> std::string_view {
    switch (value) {
    case ExecutionMode::full:
        return "full";
    case ExecutionMode::cutoff:
        return "cutoff";
    case ExecutionMode::cover:
        return "cover";
    }
    throw std::invalid_argument{"unknown execution mode"};
}

ExecutionPolicy::ExecutionPolicy(const ExecutionMode mode, const std::optional<double> radius)
    : mode_{mode}, radius_{radius} {
    validate_policy(mode_, radius_);
}

auto ExecutionPolicy::mode() const noexcept -> ExecutionMode {
    return mode_;
}

auto ExecutionPolicy::radius() const noexcept -> std::optional<double> {
    return radius_;
}

ExecutionSelection::ExecutionSelection(const ExecutionSelectionKind kind,
                                       const std::optional<double> radius)
    : kind_{kind}, radius_{radius} {
    validate_selection(kind_, radius_);
}

auto ExecutionSelection::kind() const noexcept -> ExecutionSelectionKind {
    return kind_;
}

auto ExecutionSelection::radius() const noexcept -> std::optional<double> {
    return radius_;
}

} // namespace chargefw::calculation
