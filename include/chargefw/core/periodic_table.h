#pragma once

#include <optional>
#include <span>
#include <string_view>

namespace chargefw::core {

struct Element {
    int atomic_number = 0;
    std::string_view symbol;
    std::string_view name;

    int period = 0;
    int group = 0;

    double covalent_radius = 0.0;
    double van_der_waals_radius = 0.0;
    double electronegativity = 0.0;
    double electron_affinity = 0.0;
    double first_ionization_potential = 0.0;

    [[nodiscard]] auto valence_electron_count() const noexcept -> std::optional<int>;
};

class PeriodicTable {
  public:
    [[nodiscard]] auto elements() const noexcept -> std::span<const Element>;

    [[nodiscard]] auto element(int atomic_number) const -> const Element&;
    [[nodiscard]] auto element(std::string_view symbol) const -> const Element&;

    [[nodiscard]] auto contains(int atomic_number) const noexcept -> bool;
    [[nodiscard]] auto contains(std::string_view symbol) const noexcept -> bool;
};

[[nodiscard]] auto periodic_table() noexcept -> const PeriodicTable&;

[[nodiscard]] auto element_symbol(int atomic_number) -> std::string_view;

} // namespace chargefw::core
