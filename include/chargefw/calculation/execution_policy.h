#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>

namespace chargefw::calculation {

inline constexpr std::size_t default_full_atom_threshold = 20'000;
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

enum class ChargeCorrectionPolicy : std::uint8_t {
    none,
    uniform,
};

[[nodiscard]] auto execution_selection_kind_from_string(std::string_view value)
    -> ExecutionSelectionKind;
[[nodiscard]] auto charge_correction_policy_from_string(std::string_view value)
    -> ChargeCorrectionPolicy;
[[nodiscard]] auto to_string(ExecutionSelectionKind value) -> std::string_view;
[[nodiscard]] auto to_string(ExecutionMode value) -> std::string_view;
[[nodiscard]] auto to_string(ChargeCorrectionPolicy value) -> std::string_view;

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
    // Zero delegates the worker count to the oneTBB runtime.
    std::size_t max_threads = 0;
};

} // namespace chargefw::calculation
