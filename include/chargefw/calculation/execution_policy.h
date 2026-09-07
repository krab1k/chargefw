#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>

namespace chargefw::calculation {

inline constexpr std::size_t default_cutoff_atom_threshold = 20'000;
inline constexpr std::size_t default_cover_atom_threshold = 80'000;
inline constexpr double minimum_reduced_radius = 8.0;
inline constexpr double default_automatic_reduced_radius = 12.0;

enum class ExecutionMode : std::uint8_t {
    full,
    cutoff,
    cover,
};

enum class ExecutionSelectionKind : std::uint8_t {
    automatic,
    full,
    cutoff,
    cover,
};

[[nodiscard]] auto execution_selection_kind_from_string(std::string_view value)
    -> ExecutionSelectionKind;
[[nodiscard]] auto to_string(ExecutionSelectionKind value) -> std::string_view;
[[nodiscard]] auto to_string(ExecutionMode value) -> std::string_view;

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
    // Automatic execution promotes expensive full calculations to cutoff above this threshold.
    // nullopt means unlimited.
    std::optional<std::size_t> cutoff_atom_threshold = default_cutoff_atom_threshold;
    // Automatic execution promotes cutoff calculations to cover above this threshold. nullopt means
    // unlimited cutoff execution.
    std::optional<std::size_t> cover_atom_threshold = default_cover_atom_threshold;
    // Zero delegates the worker count to the oneTBB runtime.
    std::size_t max_threads = 0;
};

} // namespace chargefw::calculation
