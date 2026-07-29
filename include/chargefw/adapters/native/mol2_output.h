#pragma once

#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule.h>

#include <iosfwd>
#include <string>

namespace chargefw::adapters::native::mol2_output {

// Copies a MOL2 source file while replacing or adding the partial-charge field for atom records in
// the selected molecule record. All unrelated source bytes are retained.
class Mol2Writer {
  public:
    explicit Mol2Writer(std::ostream& output);

    auto write_preserving_source(const std::string& source_path, std::size_t record_index,
                                 const charges::ChargeAssignment& assignment) const -> void;

    // Generates a MOL2 record from native graph and conformer data. Generated atom types are
    // element symbols only; no Tripos typing or substructure inference is performed.
    auto write_generated(const core::Molecule& molecule,
                         const charges::ChargeAssignment& assignment) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::native::mol2_output
