#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/topology_features.h>

namespace chargefw::features {

auto PreparedMolecule::molecule() const noexcept -> const core::Molecule& {
    return topology_.molecule();
}

auto PreparedMolecule::topology() const noexcept -> const TopologyFeatures& {
    return topology_;
}

} // namespace chargefw::features