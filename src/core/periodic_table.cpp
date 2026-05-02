#include <chargefw/core/periodic_table.h>

#include "core/elements.h"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace chargefw::core {

auto Element::valence_electron_count() const noexcept -> std::optional<int> {
    if (atomic_number == 2) {
        return 2;
    }

    if (group == 1 || group == 2) {
        return group;
    }

    if (group >= 13 && group <= 18) {
        return group - 10;
    }

    return std::nullopt;
}

auto PeriodicTable::elements() const noexcept -> std::span<const Element> {
    return detail::elements();
}

auto PeriodicTable::element(const int atomic_number) const -> const Element& {
    const auto all_elements = elements();

    if (atomic_number <= 0 || atomic_number > static_cast<int>(all_elements.size())) {
        throw std::out_of_range{"atomic number is outside the bundled periodic table"};
    }

    const auto& found = all_elements[static_cast<std::size_t>(atomic_number - 1)];

    if (found.atomic_number != atomic_number) {
        throw std::logic_error{"periodic table data are not indexed by atomic number"};
    }

    return found;
}

auto PeriodicTable::element(const std::string_view symbol) const -> const Element& {
    const auto all_elements = elements();

    const auto found = std::ranges::find_if(all_elements, [symbol](const Element& element) -> bool {
        return element.symbol == symbol;
    });

    if (found == all_elements.end()) {
        throw std::out_of_range{"unknown element symbol '" + std::string{symbol} + "'"};
    }

    return *found;
}

auto PeriodicTable::contains(const int atomic_number) const noexcept -> bool {
    const auto all_elements = elements();

    return atomic_number > 0 && atomic_number <= static_cast<int>(all_elements.size()) &&
           all_elements[static_cast<std::size_t>(atomic_number - 1)].atomic_number == atomic_number;
}

auto PeriodicTable::contains(const std::string_view symbol) const noexcept -> bool {
    const auto all_elements = elements();

    return std::ranges::any_of(all_elements, [symbol](const Element& element) -> bool {
        return element.symbol == symbol;
    });
}

auto periodic_table() noexcept -> const PeriodicTable& {
    static constexpr PeriodicTable table;
    return table;
}

auto element_symbol(const int atomic_number) -> std::string_view {
    return periodic_table().element(atomic_number).symbol;
}

} // namespace chargefw::core
