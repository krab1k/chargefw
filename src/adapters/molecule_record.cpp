#include <chargefw/adapters/molecule_record.h>

#include <cstddef>

namespace chargefw::adapters {

auto is_identity_mapping(const MoleculeRecordMapping& mapping,
                         const core::Molecule& molecule) noexcept -> bool {
    if (mapping.atom_indices.size() != molecule.atom_count() ||
        mapping.conformer_indices.size() != molecule.conformer_count()) {
        return false;
    }

    for (std::size_t index = 0; index < mapping.atom_indices.size(); ++index) {
        if (mapping.atom_indices[index] != index) {
            return false;
        }
    }

    for (std::size_t index = 0; index < mapping.conformer_indices.size(); ++index) {
        if (mapping.conformer_indices[index] != index) {
            return false;
        }
    }

    return true;
}

} // namespace chargefw::adapters
