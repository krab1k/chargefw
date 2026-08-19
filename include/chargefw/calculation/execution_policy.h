#pragma once

#include <cstddef>
#include <optional>

namespace chargefw::calculation {

inline constexpr std::size_t default_full_atom_threshold = 20'000;
inline constexpr double minimum_reduced_radius = 8.0;

enum class ExecutionMode {
    full,
    cutoff,
    cover,
};

enum class ExecutionSelectionKind {
    automatic,
    full,
    cutoff,
    cover,
};

enum class ChargeCorrectionPolicy {
    none,
    uniform,
};

class ExecutionPolicy {
  public:
    explicit ExecutionPolicy(
        ExecutionMode mode = ExecutionMode::full, std::optional<double> radius = {},
        ChargeCorrectionPolicy charge_correction = ChargeCorrectionPolicy::none);

    [[nodiscard]] auto mode() const noexcept -> ExecutionMode;
    [[nodiscard]] auto radius() const noexcept -> std::optional<double>;
    [[nodiscard]] auto charge_correction() const noexcept -> ChargeCorrectionPolicy;

  private:
    ExecutionMode mode_ = ExecutionMode::full;
    std::optional<double> radius_;
    ChargeCorrectionPolicy charge_correction_ = ChargeCorrectionPolicy::none;
};

class ExecutionSelection {
  public:
    explicit ExecutionSelection(ExecutionSelectionKind kind = ExecutionSelectionKind::automatic,
                                std::optional<double> radius = {},
                                std::optional<ChargeCorrectionPolicy> charge_correction = {});

    [[nodiscard]] auto kind() const noexcept -> ExecutionSelectionKind;
    [[nodiscard]] auto radius() const noexcept -> std::optional<double>;
    [[nodiscard]] auto charge_correction() const noexcept -> std::optional<ChargeCorrectionPolicy>;

  private:
    ExecutionSelectionKind kind_ = ExecutionSelectionKind::automatic;
    std::optional<double> radius_;
    std::optional<ChargeCorrectionPolicy> charge_correction_;
};

struct ResourcePolicy {
    // nullopt means unlimited.
    std::optional<std::size_t> full_atom_threshold = default_full_atom_threshold;
};

} // namespace chargefw::calculation
