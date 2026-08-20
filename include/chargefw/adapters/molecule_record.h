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

struct MoleculeRecordDiagnostic {
    std::string message;
    std::optional<std::size_t> line;
};

// A successful import result. Molecule ownership and source identity travel together so output
// adapters can retain source ordering and identity.
struct ImportedMoleculeRecord {
    core::Molecule molecule;
    MoleculeRecordIdentity identity;
    std::vector<MoleculeRecordDiagnostic> diagnostics;
};

} // namespace chargefw::adapters
