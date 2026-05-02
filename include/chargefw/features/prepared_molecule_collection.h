#pragma once

#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule.h>

#include <cstddef>
#include <span>
#include <vector>

namespace chargefw::features {

class PreparedMoleculeCollection {
  public:
    explicit PreparedMoleculeCollection(const core::MoleculeCollection& collection) {
        molecules_.reserve(collection.size());

        for (const auto& molecule : collection.molecules()) {
            molecules_.emplace_back(molecule);
        }
    }

    PreparedMoleculeCollection(core::MoleculeCollection&&) = delete;
    PreparedMoleculeCollection(const core::MoleculeCollection&&) = delete;

    [[nodiscard]] auto molecules() const noexcept -> std::span<const PreparedMolecule> {
        return molecules_;
    }

    [[nodiscard]] auto molecule_count() const noexcept -> std::size_t {
        return molecules_.size();
    }

    [[nodiscard]] auto empty() const noexcept -> bool {
        return molecules_.empty();
    }

    [[nodiscard]] auto operator[](std::size_t index) const noexcept -> const PreparedMolecule& {
        return molecules_[index];
    }

    [[nodiscard]] auto at(std::size_t index) const -> const PreparedMolecule& {
        return molecules_.at(index);
    }

  private:
    std::vector<PreparedMolecule> molecules_;
};

} // namespace chargefw::features