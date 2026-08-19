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

class ExecutionPolicy {
  public:
    explicit ExecutionPolicy(ExecutionMode mode = ExecutionMode::full,
                             std::optional<double> radius = {});

    [[nodiscard]] auto mode() const noexcept -> ExecutionMode;
    [[nodiscard]] auto radius() const noexcept -> std::optional<double>;

  private:
    ExecutionMode mode_ = ExecutionMode::full;
    std::optional<double> radius_;
};

class ExecutionSelection {
  public:
    explicit ExecutionSelection(ExecutionSelectionKind kind = ExecutionSelectionKind::automatic,
                                std::optional<double> radius = {});

    [[nodiscard]] auto kind() const noexcept -> ExecutionSelectionKind;
    [[nodiscard]] auto radius() const noexcept -> std::optional<double>;

  private:
    ExecutionSelectionKind kind_ = ExecutionSelectionKind::automatic;
    std::optional<double> radius_;
};

struct ResourcePolicy {
    // nullopt means unlimited.
    std::optional<std::size_t> full_atom_threshold = default_full_atom_threshold;
};

} // namespace chargefw::calculation
