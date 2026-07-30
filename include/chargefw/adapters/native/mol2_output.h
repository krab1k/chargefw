#pragma once

#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule.h>

#include <iosfwd>
#include <span>
#include <string>

namespace chargefw::adapters::native::mol2_output {

// Copies a MOL2 source record while replacing or adding partial-charge fields. All unrelated source
// bytes are retained.
class Mol2Writer {
  public:
    explicit Mol2Writer(std::ostream& output);

    // Streams a MOL2 file once, preserving every record and patching assignments by molecule index.
    auto write_preserving_source(const std::string& source_path,
                                 std::span<const charges::ChargeAssignment> assignments) const
        -> void;

    // Preserves MOL2 source bytes while patching assignments by molecule index.
    auto write_preserving_buffer(std::string_view source,
                                 std::span<const charges::ChargeAssignment> assignments) const
        -> void;

    // Generates a MOL2 record from native graph and conformer data. Generated atom types are
    // element symbols only; no Tripos typing or substructure inference is performed.
    auto write_generated(const core::Molecule& molecule,
                         const charges::ChargeAssignment& assignment) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::native::mol2_output
