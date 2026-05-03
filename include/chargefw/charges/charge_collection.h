#pragma once

#include <chargefw/charges/atomic_charges.h>

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::charges {

struct ChargeTarget {
    std::size_t molecule_index = 0;
    std::optional<std::size_t> conformer_index;

    [[nodiscard]] auto is_conformer_specific() const noexcept -> bool {
        return conformer_index.has_value();
    }
};

struct ChargeAssignment {
    ChargeTarget target;
    AtomicCharges charges;
};

class ChargeSet {
  public:
    explicit ChargeSet(std::string method_id, std::vector<ChargeAssignment> assignments,
                       std::optional<std::string> parameter_set_id = std::nullopt);

    [[nodiscard]] auto method_id() const noexcept -> std::string_view;

    [[nodiscard]] auto parameter_set_id() const noexcept -> std::optional<std::string_view>;

    [[nodiscard]] auto assignments() const noexcept -> std::span<const ChargeAssignment>;

    [[nodiscard]] auto assignment(std::size_t index) const -> const ChargeAssignment&;

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto empty() const noexcept -> bool;

  private:
    std::string method_id_;
    std::optional<std::string> parameter_set_id_;
    std::vector<ChargeAssignment> assignments_;
};

class ChargeCollection {
  public:
    explicit ChargeCollection(std::vector<ChargeSet> charge_sets);

    [[nodiscard]] auto charge_sets() const noexcept -> std::span<const ChargeSet>;

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto empty() const noexcept -> bool;

    [[nodiscard]] auto operator[](std::size_t index) const noexcept -> const ChargeSet&;
    [[nodiscard]] auto at(std::size_t index) const -> const ChargeSet&;

  private:
    std::vector<ChargeSet> charge_sets_;
};

} // namespace chargefw::charges