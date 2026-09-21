#pragma once

#include <chargefw/core/molecule.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::adapters {

enum class MolecularSourceFormat : std::uint8_t {
    molecule_json,
    mol,
    sdf,
    mol2,
    pdb,
    mmcif,
};

enum class SourceConnectivity : std::uint8_t {
    absent,
    explicitly_empty,
    present,
};

struct SourceHierarchyLabels {
    std::optional<std::string> atom;
    std::optional<std::string> residue;
    std::optional<std::string> chain;
    std::optional<std::string> sequence;

    [[nodiscard]] auto operator==(const SourceHierarchyLabels&) const -> bool = default;
};

struct SourceStructuralLabels {
    SourceHierarchyLabels author;
    SourceHierarchyLabels label;
    std::optional<std::string> entity;
    std::optional<std::string> insertion_code;
    std::optional<std::string> alternate_location;
    std::optional<std::string> segment;

    [[nodiscard]] auto operator==(const SourceStructuralLabels&) const -> bool = default;
};

// Source positions are zero-based within one record. IDs retain explicit source tokens when the
// format provides them; position remains authoritative when no portable ID exists.
struct SourceAtomReference {
    std::size_t position = 0;
    std::optional<std::string> id;
    std::optional<SourceStructuralLabels> structural_labels;

    [[nodiscard]] auto operator==(const SourceAtomReference&) const -> bool = default;
};

struct SourceConformerReference {
    std::size_t position = 0;
    std::optional<std::string> id;
    std::vector<SourceAtomReference> sites;

    [[nodiscard]] auto operator==(const SourceConformerReference&) const -> bool = default;
};

// Small adapter-owned metadata retained with the normalized molecule. Structural labels are added
// by the structural-reader boundary rather than by core::Molecule.
struct MoleculeImportMetadata {
    MolecularSourceFormat format = MolecularSourceFormat::molecule_json;
    std::vector<SourceAtomReference> atoms;
    std::vector<SourceConformerReference> conformers;
    std::optional<std::string> record_selection;
    std::optional<std::string> alternate_location_selection;
    std::optional<std::string> conformer_selection;
    std::optional<std::string> bond_strategy;
    SourceConnectivity source_connectivity = SourceConnectivity::absent;
};

// Identifies one source record without imposing a file-format or toolkit dependency.
struct MoleculeRecordIdentity {
    std::string source;
    std::size_t record_index = 0;
    std::string record_id;
};

struct MoleculeRecordDiagnostic {
    std::string code;
    std::string message;
    std::optional<std::size_t> line;
};

// A successful import result. Molecule ownership and source identity travel together so output
// adapters can retain source ordering and identity.
struct ImportedMoleculeRecord {
    core::Molecule molecule;
    MoleculeRecordIdentity identity;
    std::vector<MoleculeRecordDiagnostic> diagnostics;
    std::optional<MoleculeImportMetadata> import_metadata;
};

} // namespace chargefw::adapters
