#pragma once

#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule.h>

#include <cstdint>
#include <iosfwd>
#include <span>
#include <string>
#include <string_view>

namespace chargefw::adapters::native::sdf_output {

enum class WriteMode : std::uint8_t { replace, append };
enum class MolFormat : std::uint8_t { v2000, v3000 };

struct ChargeProperty {
    std::size_t charge_type_id = 0;
    std::span<const charges::ChargeAssignment> assignments;
    std::string_view method;
    std::string_view parameter_set;
    std::string_view software_name;
    std::string_view software_version;
};

// Copies an SDF source while attaching ChargeFW charge properties and their method metadata before
// each record delimiter. Replace mode removes existing CHARGEFW_CHARGES_* and
// CHARGEFW_CHARGE_METADATA_* fields first; append mode retains them.
class SdfWriter {
  public:
    explicit SdfWriter(std::ostream& output);

    auto write_preserving_source(const std::string& source_path,
                                 std::span<const ChargeProperty> properties,
                                 WriteMode mode = WriteMode::replace) const -> void;

    auto write_preserving_buffer(std::string_view source,
                                 std::span<const ChargeProperty> properties,
                                 WriteMode mode = WriteMode::replace) const -> void;

    // Generates one SDF record from native graph and conformer data. Partial charges are written as
    // numbered properties; formal charges remain part of the selected MOL representation.
    auto write_generated(const core::Molecule& molecule, std::span<const ChargeProperty> properties,
                         MolFormat format = MolFormat::v2000) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::native::sdf_output
