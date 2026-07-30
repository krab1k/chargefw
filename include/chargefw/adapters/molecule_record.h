#pragma once

#include <chargefw/core/molecule.h>

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::adapters {

// Identifies one source record without imposing a file-format or toolkit dependency.
struct MoleculeRecordIdentity {
    std::string source;
    std::size_t record_index = 0;
    std::string record_id;
};

// Maps a source atom or conformer index to its preserved native index. A missing value represents
// an explicitly omitted source item; adapters must not use it to silently reorder chemistry.
struct MoleculeRecordMapping {
    std::vector<std::optional<std::size_t>> atom_indices;
    std::vector<std::optional<std::size_t>> conformer_indices;
};

struct MoleculeRecordDiagnostic {
    std::string message;
    std::optional<std::size_t> line;
};

// A successful import result. Molecule ownership, source identity, and source-to-native mapping
// travel together so output adapters can retain source ordering and identity.
struct ImportedMoleculeRecord {
    core::Molecule molecule;
    MoleculeRecordIdentity identity;
    MoleculeRecordMapping mapping;
    std::vector<MoleculeRecordDiagnostic> diagnostics;
};

[[nodiscard]] auto is_identity_mapping(const MoleculeRecordMapping& mapping,
                                       const core::Molecule& molecule) noexcept -> bool;

} // namespace chargefw::adapters
