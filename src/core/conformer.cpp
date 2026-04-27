#include <chargefw/core/conformer.h>

#include <utility>

namespace chargefw::core {

Conformer::Conformer(std::vector<Position> positions, std::string name)
    : positions_{std::move(positions)}, name_{std::move(name)} {}

auto Conformer::positions() const noexcept -> std::span<const Position> {
    return positions_;
}

auto Conformer::position(const std::size_t atom_index) const -> const Position& {
    return positions_.at(atom_index);
}

auto Conformer::size() const noexcept -> std::size_t {
    return positions_.size();
}

auto Conformer::empty() const noexcept -> bool {
    return positions_.empty();
}

auto Conformer::name() const noexcept -> std::string_view {
    return name_;
}

auto Conformer::operator[](std::size_t index) const noexcept -> const Position& {
        return positions_[index];
}

auto Conformer::at(std::size_t index) const -> const Position& {
    return positions_.at(index);
}

} // namespace chargefw::core