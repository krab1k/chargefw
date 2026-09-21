#pragma once

#include <chargefw/adapters/molecule_record.h>

#include <optional>
#include <utility>
#include <vector>

namespace chargefw::python {

class NativeInputMetadata {
  public:
    NativeInputMetadata(adapters::MoleculeRecordIdentity identity,
                        std::vector<adapters::MoleculeRecordDiagnostic> diagnostics,
                        std::optional<adapters::MoleculeImportMetadata> import_metadata)
        : identity_{std::move(identity)}, diagnostics_{std::move(diagnostics)},
          import_metadata_{std::move(import_metadata)} {}

    [[nodiscard]] auto make_record(core::Molecule molecule) const
        -> adapters::ImportedMoleculeRecord {
        return adapters::ImportedMoleculeRecord{.molecule = std::move(molecule),
                                                .identity = identity_,
                                                .diagnostics = diagnostics_,
                                                .import_metadata = import_metadata_};
    }

  private:
    adapters::MoleculeRecordIdentity identity_;
    std::vector<adapters::MoleculeRecordDiagnostic> diagnostics_;
    std::optional<adapters::MoleculeImportMetadata> import_metadata_;
};

} // namespace chargefw::python
