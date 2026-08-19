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

auto validate_charge_correction(const ChargeCorrectionPolicy charge_correction) -> void {
    switch (charge_correction) {
    case ChargeCorrectionPolicy::none:
    case ChargeCorrectionPolicy::uniform:
        return;
    }

    throw std::invalid_argument{"unknown charge correction policy"};
}

auto validate_policy(const ExecutionMode mode, const std::optional<double> radius,
                     const ChargeCorrectionPolicy charge_correction) -> void {
    validate_charge_correction(charge_correction);

    switch (mode) {
    case ExecutionMode::full:
        if (radius.has_value()) {
            throw std::invalid_argument{"full execution does not accept a radius"};
        }
        if (charge_correction != ChargeCorrectionPolicy::none) {
            throw std::invalid_argument{"full execution does not accept a charge correction"};
        }
        return;
    case ExecutionMode::cutoff:
    case ExecutionMode::cover:
        validate_reduced_radius(radius);
        return;
    }

    throw std::invalid_argument{"unknown execution mode"};
}

auto validate_selection(const ExecutionSelectionKind kind, const std::optional<double> radius,
                        const std::optional<ChargeCorrectionPolicy> charge_correction) -> void {
    if (charge_correction.has_value()) {
        validate_charge_correction(*charge_correction);
    }

    switch (kind) {
    case ExecutionSelectionKind::automatic:
        if (radius.has_value()) {
            validate_reduced_radius(radius);
        }
        if (charge_correction.has_value()) {
            throw std::invalid_argument{"automatic execution does not accept a charge correction"};
        }
        return;
    case ExecutionSelectionKind::full:
        validate_policy(ExecutionMode::full, radius, ChargeCorrectionPolicy::none);
        if (charge_correction.has_value()) {
            throw std::invalid_argument{"full execution does not accept a charge correction"};
        }
        return;
    case ExecutionSelectionKind::cutoff:
        validate_policy(ExecutionMode::cutoff, radius,
                        charge_correction.value_or(ChargeCorrectionPolicy::uniform));
        return;
    case ExecutionSelectionKind::cover:
        validate_policy(ExecutionMode::cover, radius,
                        charge_correction.value_or(ChargeCorrectionPolicy::uniform));
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

auto charge_correction_policy_from_string(const std::string_view value) -> ChargeCorrectionPolicy {
    if (value == "none") {
        return ChargeCorrectionPolicy::none;
    }
    if (value == "uniform") {
        return ChargeCorrectionPolicy::uniform;
    }
    throw std::invalid_argument{"unknown charge correction policy: " + std::string{value}};
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

auto to_string(const ChargeCorrectionPolicy value) -> std::string_view {
    switch (value) {
    case ChargeCorrectionPolicy::none:
        return "none";
    case ChargeCorrectionPolicy::uniform:
        return "uniform";
    }
    throw std::invalid_argument{"unknown charge correction policy"};
}

ExecutionPolicy::ExecutionPolicy(const ExecutionMode mode, const std::optional<double> radius,
                                 const ChargeCorrectionPolicy charge_correction)
    : mode_{mode}, radius_{radius}, charge_correction_{charge_correction} {
    validate_policy(mode_, radius_, charge_correction_);
}

auto ExecutionPolicy::mode() const noexcept -> ExecutionMode {
    return mode_;
}

auto ExecutionPolicy::radius() const noexcept -> std::optional<double> {
    return radius_;
}

auto ExecutionPolicy::charge_correction() const noexcept -> ChargeCorrectionPolicy {
    return charge_correction_;
}

ExecutionSelection::ExecutionSelection(
    const ExecutionSelectionKind kind, const std::optional<double> radius,
    const std::optional<ChargeCorrectionPolicy> charge_correction)
    : kind_{kind}, radius_{radius}, charge_correction_{charge_correction} {
    validate_selection(kind_, radius_, charge_correction_);
}

auto ExecutionSelection::kind() const noexcept -> ExecutionSelectionKind {
    return kind_;
}

auto ExecutionSelection::radius() const noexcept -> std::optional<double> {
    return radius_;
}

auto ExecutionSelection::charge_correction() const noexcept
    -> std::optional<ChargeCorrectionPolicy> {
    return charge_correction_;
}

} // namespace chargefw::calculation
