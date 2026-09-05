#include <chargefw/adapters/generated_output.h>

#include <chargefw/adapters/gemmi/mmcif_output.h>
#include <chargefw/adapters/native/mol2_output.h>
#include <chargefw/adapters/native/sdf_output.h>

#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::adapters::generated_output {
namespace {

auto validate_conformer(const core::Molecule& molecule, const std::size_t conformer_index) -> void {
    if (conformer_index >= molecule.conformer_count()) {
        throw std::invalid_argument{"selected output conformer is unavailable"};
    }
    for (const auto& position : molecule.conformer(conformer_index).positions()) {
        if (!std::isfinite(position.x) || !std::isfinite(position.y) ||
            !std::isfinite(position.z)) {
            throw std::invalid_argument{"molecular output coordinates must be finite"};
        }
    }
}

auto selected_charge_set(const charges::ChargeSet& charge_set,
                         const std::span<const ImportedMoleculeRecord> records)
    -> charges::ChargeSet {
    auto selected = std::vector<charges::ChargeAssignment>{};
    selected.reserve(records.size());
    for (std::size_t molecule_index = 0; molecule_index < records.size(); ++molecule_index) {
        const charges::ChargeAssignment* selected_assignment = nullptr;
        for (const auto& candidate : charge_set.assignments()) {
            if (candidate.target.molecule_index == molecule_index &&
                (!candidate.target.conformer_index.has_value() ||
                 candidate.target.conformer_index == 0)) {
                selected_assignment = std::addressof(candidate);
                if (selected_assignment->target.conformer_index == 0) {
                    break;
                }
            }
        }
        if (selected_assignment == nullptr) {
            throw std::invalid_argument{
                "first conformer has no charge assignment for output molecule " +
                std::to_string(molecule_index)};
        }
        const auto& molecule = records[molecule_index].molecule;
        validate_conformer(molecule, 0);
        selected.push_back(charges::ChargeAssignment{
            .target = {.molecule_index = molecule_index, .conformer_index = 0},
            .charges = selected_assignment->charges});
    }
    return charges::ChargeSet{std::string{charge_set.method_id()}, std::move(selected),
                              charge_set.parameter_set_id().transform(
                                  [](const std::string_view id) { return std::string{id}; })};
}

auto write_mmcif(std::ostream& output, const std::span<const ImportedMoleculeRecord> records,
                 const charges::ChargeSet& charge_set, const std::string_view generator_name,
                 const std::string_view generator_version) -> void {
    for (const auto& record : records) {
        if (record.molecule.conformer_count() == 0) {
            throw std::invalid_argument{"mmCIF output requires coordinates"};
        }
        for (std::size_t conformer_index = 0; conformer_index < record.molecule.conformer_count();
             ++conformer_index) {
            validate_conformer(record.molecule, conformer_index);
        }
    }
    gemmi::mmcif_output::MmcifWriter{output}.write_generated(records, charge_set, generator_name,
                                                             generator_version);
}

auto write_single_conformer(std::ostream& output,
                            const std::span<const ImportedMoleculeRecord> records,
                            const charges::ChargeSet& charge_set, const Format format,
                            const std::string_view generator_name,
                            const std::string_view generator_version) -> void {
    const auto selected = selected_charge_set(charge_set, records);
    for (std::size_t index = 0; index < records.size(); ++index) {
        const auto assignment = selected.assignments().subspan(index, 1);
        if (format == Format::mol2) {
            native::mol2_output::Mol2Writer{output}.write_generated(records[index].molecule,
                                                                    assignment.front());
            continue;
        }
        const auto property = std::array{native::sdf_output::ChargeProperty{
            .charge_type_id = 1,
            .assignments = assignment,
            .method = selected.method_id(),
            .parameter_set = selected.parameter_set_id().value_or(""),
            .software_name = generator_name,
            .software_version = generator_version}};
        native::sdf_output::SdfWriter{output}.write_generated(
            records[index].molecule, property,
            format == Format::sdf_v2000 ? native::sdf_output::MolFormat::v2000
                                        : native::sdf_output::MolFormat::v3000);
    }
}

} // namespace

auto write(std::ostream& output, const std::span<const ImportedMoleculeRecord> records,
           const charges::ChargeSet& charge_set, const Format format,
           const std::string_view generator_name, const std::string_view generator_version)
    -> void {
    switch (format) {
    case Format::sdf_v2000:
    case Format::sdf_v3000:
    case Format::mol2:
        write_single_conformer(output, records, charge_set, format, generator_name,
                               generator_version);
        return;
    case Format::mmcif:
        write_mmcif(output, records, charge_set, generator_name, generator_version);
        return;
    }
    throw std::invalid_argument{"unknown generated molecular output format"};
}

} // namespace chargefw::adapters::generated_output
