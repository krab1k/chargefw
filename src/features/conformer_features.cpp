#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>

#include <cmath>
#include <stdexcept>

namespace chargefw::features {
namespace {

auto validate_radius(const double radius) -> void {
    if (radius < 0.0) {
        throw std::invalid_argument{"radius must not be negative"};
    }
}

} // namespace

ConformerFeatures::ConformerFeatures(const core::Molecule& molecule,
                                     const std::size_t conformer_index)
    : molecule_{&molecule}, conformer_index_{conformer_index} {
    validate_conformer_index();
}

auto ConformerFeatures::molecule() const noexcept -> const core::Molecule& {
    return *molecule_;
}

auto ConformerFeatures::conformer_index() const noexcept -> std::size_t {
    return conformer_index_;
}

auto ConformerFeatures::position(const std::size_t atom_index) const -> const core::Position& {
    validate_atom_index(atom_index);
    return molecule_->conformer(conformer_index_)[atom_index];
}

auto ConformerFeatures::squared_distance(const std::size_t first_atom_index,
                                         const std::size_t second_atom_index) const -> double {
    return core::squared_distance(position(first_atom_index), position(second_atom_index));
}

auto ConformerFeatures::distance(const std::size_t first_atom_index,
                                 const std::size_t second_atom_index) const -> double {
    return core::distance(position(first_atom_index), position(second_atom_index));
}

auto ConformerFeatures::neighbor_indices_within(const std::size_t atom_index,
                                                const double radius) const
    -> std::vector<std::size_t> {
    validate_atom_index(atom_index);
    validate_radius(radius);

    const auto radius_squared = radius * radius;
    std::vector<std::size_t> neighbors;

    for (std::size_t other_atom_index = 0; other_atom_index < molecule_->atom_count();
         ++other_atom_index) {
        if (other_atom_index == atom_index) {
            continue;
        }

        if (squared_distance(atom_index, other_atom_index) <= radius_squared) {
            neighbors.push_back(other_atom_index);
        }
    }

    return neighbors;
}

auto ConformerFeatures::validate_conformer_index() const -> void {
    if (conformer_index_ >= molecule_->conformer_count()) {
        throw std::out_of_range{"conformer index is outside the molecule"};
    }
}

auto ConformerFeatures::validate_atom_index(const std::size_t atom_index) const -> void {
    if (atom_index >= molecule_->atom_count()) {
        throw std::out_of_range{"atom index is outside the molecule"};
    }
}

} // namespace chargefw::features
