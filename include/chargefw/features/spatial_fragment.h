#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/parameters/classification/parameter_classification.h>

#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace chargefw::features {

inline constexpr std::size_t no_source_index = std::numeric_limits<std::size_t>::max();

class SpatialFragment {
  public:
    [[nodiscard]] auto molecule() const noexcept -> const core::Molecule&;
    [[nodiscard]] auto source_conformer_index() const noexcept -> std::size_t;
    [[nodiscard]] auto center_local_atom_index() const noexcept -> std::size_t;

    // Local atom order follows the source molecule's atom order.
    [[nodiscard]] auto local_to_source_atom_indices() const noexcept
        -> std::span<const std::size_t>;
    // Unselected source atoms have no_source_index. Reassess this source-sized lookup before
    // retaining many fragments or adding parallel fragment execution.
    [[nodiscard]] auto source_to_local_atom_indices() const noexcept
        -> std::span<const std::size_t>;
    [[nodiscard]] auto local_to_source_bond_indices() const noexcept
        -> std::span<const std::size_t>;

  private:
    SpatialFragment(core::Molecule molecule, std::vector<std::size_t> local_to_source_atom_indices,
                    std::vector<std::size_t> source_to_local_atom_indices,
                    std::vector<std::size_t> local_to_source_bond_indices,
                    std::size_t center_local_atom_index, std::size_t source_conformer_index);

    core::Molecule molecule_;
    std::vector<std::size_t> local_to_source_atom_indices_;
    std::vector<std::size_t> source_to_local_atom_indices_;
    std::vector<std::size_t> local_to_source_bond_indices_;
    std::size_t center_local_atom_index_;
    std::size_t source_conformer_index_;

    friend auto build_spatial_fragment(const PreparedMolecule&, std::size_t, std::size_t, double)
        -> SpatialFragment;
};

[[nodiscard]] auto build_spatial_fragment(const PreparedMolecule& source,
                                          std::size_t conformer_index,
                                          std::size_t center_atom_index, double radius)
    -> SpatialFragment;

// Projects whole-molecule parameter entry indices to an induced fragment without reclassification.
[[nodiscard]] auto project_classification(const parameters::ParameterClassification& source,
                                          const SpatialFragment& fragment)
    -> parameters::ParameterClassification;

} // namespace chargefw::features
