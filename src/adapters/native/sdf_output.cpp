#include <chargefw/adapters/native/sdf_output.h>

#include "common_output.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <format>
#include <fstream>
#include <optional>
#include <ostream>
#include <print>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

namespace chargefw::adapters::native::sdf_output {
namespace {

constexpr std::string_view delimiter{"$$$$"};
constexpr std::string_view property_prefix{"CHARGEFW_CHARGES_"};

struct Line {
    std::string_view content;
    std::string_view ending;
    std::size_t end = 0;
};

[[nodiscard]] auto formatted_charge(const double value) -> std::string {
    constexpr auto charge_scale = 10000.0;
    return std::format("{:.4f}", std::round(value * charge_scale) / charge_scale);
}

auto validate_property_ids(const std::span<const ChargeProperty> properties) -> void {
    auto ids = std::unordered_set<std::size_t>{};
    for (const auto& property : properties) {
        if (property.charge_type_id == 0 || !ids.insert(property.charge_type_id).second) {
            throw std::invalid_argument{"SDF charge type IDs must be unique positive integers"};
        }
    }
}

[[nodiscard]] auto line_at(const std::string_view value, const std::size_t start) -> Line {
    const auto newline = value.find('\n', start);
    const auto end = newline == std::string_view::npos ? value.size() : newline + 1;
    auto content = value.substr(start, end - start);
    auto ending = std::string_view{};
    if (content.ends_with("\r\n")) {
        content.remove_suffix(2);
        ending = "\r\n";
    } else if (content.ends_with('\n')) {
        content.remove_suffix(1);
        ending = "\n";
    }
    return {.content = content, .ending = ending, .end = end};
}

[[nodiscard]] auto property_name(const std::string_view line) -> std::optional<std::string_view> {
    if (!line.starts_with('>')) {
        return std::nullopt;
    }
    const auto open = line.find('<');
    const auto close =
        open == std::string_view::npos ? std::string_view::npos : line.find('>', open);
    if (open == std::string_view::npos || close == std::string_view::npos) {
        return std::nullopt;
    }
    return line.substr(open + 1, close - open - 1);
}

[[nodiscard]] auto is_chargefw_property(const std::string_view line) -> bool {
    const auto name = property_name(line);
    if (!name.has_value() || !name->starts_with(property_prefix)) {
        return false;
    }
    const auto id = name->substr(property_prefix.size());
    if (id.empty()) {
        return false;
    }
    std::size_t parsed = 0;
    const auto [pointer, error] = std::from_chars(id.data(), id.data() + id.size(), parsed);
    return error == std::errc{} && pointer == id.data() + id.size() && parsed > 0;
}

[[nodiscard]] auto record_ending(const std::string_view record) -> std::string_view {
    std::size_t start = 0;
    while (start < record.size()) {
        const auto line = line_at(record, start);
        if (!line.ending.empty()) {
            return line.ending;
        }
        start = line.end;
    }
    return "\n";
}

auto write_source_fields(const std::string_view record, const WriteMode mode, std::ostream& output)
    -> void {
    if (mode == WriteMode::append) {
        std::print(output, "{}", record);
        return;
    }

    std::size_t start = 0;
    while (start < record.size()) {
        const auto line = line_at(record, start);
        if (!is_chargefw_property(line.content)) {
            std::print(output, "{}", record.substr(start, line.end - start));
            start = line.end;
            continue;
        }

        start = line.end;
        while (start < record.size()) {
            const auto value_line = line_at(record, start);
            start = value_line.end;
            if (value_line.content.empty()) {
                break;
            }
        }
    }
}

auto write_charge_property(const ChargeProperty& property,
                           const charges::ChargeAssignment& assignment,
                           const std::string_view ending, std::ostream& output) -> void {
    if (property.charge_type_id == 0) {
        throw std::invalid_argument{"SDF charge type IDs must start at 1"};
    }
    std::print(output, "> <{}{}>{}", property_prefix, property.charge_type_id, ending);
    for (std::size_t index = 0; index < assignment.charges.size(); ++index) {
        if (index != 0) {
            std::print(output, " ");
        }
        std::print(output, "{}", formatted_charge(assignment.charges[index]));
    }
    std::print(output, "{}{}", ending, ending);
}

auto write_record(const std::string_view record, const std::string_view ending,
                  const std::size_t record_index, const std::span<const ChargeProperty> properties,
                  const WriteMode mode, std::ostream& output) -> void {
    const auto selected_ending = ending.empty() ? record_ending(record) : ending;
    write_source_fields(record, mode, output);
    if (!record.empty() && !record.ends_with('\n')) {
        std::print(output, "{}", selected_ending);
    }

    for (const auto& property : properties) {
        if (record_index >= property.assignments.size()) {
            throw std::invalid_argument{"SDF source has more records than charge assignments"};
        }
        const auto& assignment = property.assignments[record_index];
        if (assignment.target.molecule_index != record_index) {
            throw std::invalid_argument{"SDF assignment order does not match source record order"};
        }
        if (properties.front().assignments[record_index].charges.size() !=
            assignment.charges.size()) {
            throw std::invalid_argument{
                "SDF charge properties for a record must have equal atom counts"};
        }
        write_charge_property(property, assignment, selected_ending, output);
    }
    std::print(output, "{}{}", delimiter, selected_ending);
}

auto write_preserving_records(std::istream& input, std::ostream& output,
                              const std::span<const ChargeProperty> properties,
                              const WriteMode mode) -> void {
    validate_property_ids(properties);

    std::size_t record_index = 0;
    auto record = std::string{};
    auto line = std::string{};
    auto record_line_ending = std::string_view{};
    while (std::getline(input, line)) {
        const auto carriage_return = line.ends_with('\r');
        if (record_line_ending.empty()) {
            record_line_ending = carriage_return ? "\r\n" : "\n";
        }
        const auto content = carriage_return ? std::string_view{line}.substr(0, line.size() - 1)
                                             : std::string_view{line};
        if (content == delimiter) {
            write_record(record, record_line_ending, record_index++, properties, mode, output);
            record.clear();
            record_line_ending = {};
        } else {
            if (carriage_return) {
                record.append(line.data(), line.size() - 1);
                record += "\r\n";
            } else {
                record += line;
                record += '\n';
            }
        }
    }

    if (!record.empty()) {
        throw std::runtime_error{"SDF source record is missing the $$$$ delimiter"};
    }
    for (const auto& property : properties) {
        if (record_index != property.assignments.size()) {
            throw std::invalid_argument{"SDF source has fewer records than charge assignments"};
        }
    }
}

[[nodiscard]] auto generated_assignments(const std::span<const ChargeProperty> properties,
                                         const core::Molecule& molecule)
    -> std::vector<const charges::ChargeAssignment*> {
    if (properties.empty()) {
        throw std::invalid_argument{"generated SDF output requires at least one charge property"};
    }
    validate_property_ids(properties);

    auto result = std::vector<const charges::ChargeAssignment*>{};
    result.reserve(properties.size());
    std::optional<std::size_t> conformer_index;
    for (const auto& property : properties) {
        if (property.assignments.size() != 1) {
            throw std::invalid_argument{
                "generated SDF output requires one assignment per charge property"};
        }
        const auto& assignment = property.assignments.front();
        common_output::validate_assignment(assignment.charges, molecule.atom_count());
        if (!assignment.target.conformer_index.has_value()) {
            throw std::invalid_argument{
                "generated SDF output requires conformer-specific assignments"};
        }
        if (conformer_index.has_value() && conformer_index != assignment.target.conformer_index) {
            throw std::invalid_argument{
                "generated SDF charge properties must reference the same conformer"};
        }
        conformer_index = assignment.target.conformer_index;
        result.push_back(std::addressof(assignment));
    }
    if (*conformer_index >= molecule.conformer_count()) {
        throw std::invalid_argument{"charge assignment references an unavailable conformer"};
    }
    return result;
}

auto write_v2000(const core::Molecule& molecule, const core::Conformer& conformer,
                 std::ostream& output) -> void {
    if (molecule.atom_count() > 999 || molecule.bond_count() > 999) {
        throw std::invalid_argument{"V2000 output supports at most 999 atoms and 999 bonds"};
    }
    std::print(output, "{}\n  ChargeFW\n\n{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
               molecule.name().empty() ? "chargefw" : molecule.name(), molecule.atom_count(),
               molecule.bond_count());
    for (std::size_t index = 0; index < molecule.atom_count(); ++index) {
        const auto& position = conformer[index];
        std::print(output, "{:>10.4f}{:>10.4f}{:>10.4f} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
                   position.x, position.y, position.z,
                   common_output::atom_element_symbol(molecule.atom(index)));
    }
    for (const auto& bond : molecule.bonds()) {
        std::print(output, "{:>3}{:>3}{:>3}  0  0  0  0\n", bond.first_atom_index() + 1,
                   bond.second_atom_index() + 1,
                   common_output::bond_type(bond.order(), BondFormat::mol));
    }

    constexpr std::size_t charges_per_line = 8;
    auto charged_atoms = std::vector<std::pair<std::size_t, int>>{};
    for (std::size_t index = 0; index < molecule.atom_count(); ++index) {
        if (molecule.atom(index).formal_charge() != 0) {
            charged_atoms.emplace_back(index + 1, molecule.atom(index).formal_charge());
        }
    }
    for (std::size_t start = 0; start < charged_atoms.size(); start += charges_per_line) {
        const auto count = std::min(charges_per_line, charged_atoms.size() - start);
        std::print(output, "M  CHG{:>3}", count);
        for (std::size_t index = 0; index < count; ++index) {
            const auto [atom_number, charge] = charged_atoms[start + index];
            std::print(output, "{:>4}{:>4}", atom_number, charge);
        }
        std::print(output, "\n");
    }
    std::print(output, "M  END\n");
}

auto write_v3000(const core::Molecule& molecule, const core::Conformer& conformer,
                 std::ostream& output) -> void {
    std::print(output,
               "{}\n  ChargeFW\n\n  0  0  0     0  0            999 V3000\n"
               "M  V30 BEGIN CTAB\nM  V30 COUNTS {} {} 0 0 0\nM  V30 BEGIN ATOM\n",
               molecule.name().empty() ? "chargefw" : molecule.name(), molecule.atom_count(),
               molecule.bond_count());
    for (std::size_t index = 0; index < molecule.atom_count(); ++index) {
        const auto& atom = molecule.atom(index);
        const auto& position = conformer[index];
        std::print(output, "M  V30 {} {} {} {} {} 0", index + 1,
                   common_output::atom_element_symbol(atom), position.x, position.y, position.z);
        if (atom.formal_charge() != 0) {
            std::print(output, " CHG={}", atom.formal_charge());
        }
        std::print(output, "\n");
    }
    std::print(output, "M  V30 END ATOM\nM  V30 BEGIN BOND\n");
    for (std::size_t index = 0; index < molecule.bond_count(); ++index) {
        const auto& bond = molecule.bond(index);
        std::print(output, "M  V30 {} {} {} {}\n", index + 1,
                   common_output::bond_type(bond.order(), BondFormat::mol),
                   bond.first_atom_index() + 1, bond.second_atom_index() + 1);
    }
    std::print(output, "M  V30 END BOND\nM  V30 END CTAB\nM  END\n");
}

} // namespace

SdfWriter::SdfWriter(std::ostream& output) : output_{std::addressof(output)} {}

auto SdfWriter::write_preserving_source(const std::string& source_path,
                                        const std::span<const ChargeProperty> properties,
                                        const WriteMode mode) const -> void {
    auto input = std::ifstream{source_path, std::ios::binary};
    if (!input) {
        throw std::runtime_error{"unable to open SDF source file: " + source_path};
    }
    write_preserving_records(input, *output_, properties, mode);
}

auto SdfWriter::write_preserving_buffer(const std::string_view source,
                                        const std::span<const ChargeProperty> properties,
                                        const WriteMode mode) const -> void {
    auto input = std::istringstream{std::string{source}};
    write_preserving_records(input, *output_, properties, mode);
}

auto SdfWriter::write_generated(const core::Molecule& molecule,
                                const std::span<const ChargeProperty> properties,
                                const MolFormat format) const -> void {
    const auto assignments = generated_assignments(properties, molecule);
    const auto conformer_index = *assignments.front()->target.conformer_index;
    const auto& conformer = molecule.conformer(conformer_index);
    if (conformer.size() != molecule.atom_count()) {
        throw std::invalid_argument{
            "SDF conformer coordinate count does not match molecule atom count"};
    }

    if (format == MolFormat::v2000) {
        write_v2000(molecule, conformer, *output_);
    } else {
        write_v3000(molecule, conformer, *output_);
    }
    for (std::size_t index = 0; index < properties.size(); ++index) {
        write_charge_property(properties[index], *assignments[index], "\n", *output_);
    }
    std::print(*output_, "$$$$\n");
}

} // namespace chargefw::adapters::native::sdf_output
