#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/features/topology_features.h>

namespace chargefw::features {

class PreparedMolecule {
  public:
    explicit PreparedMolecule(const core::Molecule& molecule) : topology_{molecule} {}

    PreparedMolecule(core::Molecule&&) = delete;
    PreparedMolecule(const core::Molecule&&) = delete;

    [[nodiscard]] auto molecule() const noexcept -> const core::Molecule&;

    [[nodiscard]] auto topology() const noexcept -> const TopologyFeatures&;

  private:
    TopologyFeatures topology_;
};

} // namespace chargefw::features