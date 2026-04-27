#pragma once

#include <chargefw/core/molecule.h>

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace chargefw::features {

class PreparedMolecule {
  public:
    explicit PreparedMolecule(const core::Molecule& molecule);

    [[nodiscard]] auto molecule() const noexcept -> const core::Molecule&;

    [[nodiscard]] auto degree(std::size_t atom_index) const -> std::size_t;

    [[nodiscard]] auto neighbor_indices(std::size_t atom_index) const
        -> std::span<const std::size_t>;

    [[nodiscard]] auto incident_bond_indices(std::size_t atom_index) const
        -> std::span<const std::size_t>;

    [[nodiscard]] auto bond_index_between(std::size_t first_atom_index,
                                          std::size_t second_atom_index) const
        -> std::optional<std::size_t>;

    [[nodiscard]] auto are_bonded(std::size_t first_atom_index, std::size_t second_atom_index) const
        -> bool;

  private:
    const core::Molecule* molecule_;

    std::vector<std::vector<std::size_t>> neighbor_indices_;
    std::vector<std::vector<std::size_t>> incident_bond_indices_;
};

} // namespace chargefw::features
