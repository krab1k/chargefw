#include <chargefw/adapters/native/mol2_output.h>

#include "common_output.h"

#include <cstddef>
#include <optional>
#include <ostream>
#include <print>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace chargefw::adapters::native::mol2_output {
namespace {

constexpr std::string_view molecule_marker{"@<TRIPOS>MOLECULE"};
constexpr std::string_view atom_marker{"@<TRIPOS>ATOM"};
constexpr std::string_view bond_marker{"@<TRIPOS>BOND"};

[[nodiscard]] auto field_range(const std::string_view value, const std::size_t field_index)
    -> std::optional<std::pair<std::size_t, std::size_t>> {
    std::size_t position = 0;
    for (std::size_t index = 0; index <= field_index; ++index) {
        position = value.find_first_not_of(" \t", position);
        if (position == std::string_view::npos) {
            return std::nullopt;
        }
        const auto end = value.find_first_of(" \t", position);
        if (index == field_index) {
            return std::pair{position, end == std::string_view::npos ? value.size() : end};
        }
        position = end;
    }
    return std::nullopt;
}

[[nodiscard]] auto patch_atom_line(const std::string_view content, const double charge,
                                   const std::string_view ending) -> std::string {
    if (!field_range(content, 5).has_value()) {
        throw std::runtime_error{"invalid MOL2 atom record while writing charges"};
    }

    auto result = std::string{content};
    const auto formatted = common_output::formatted_charge(charge);
    if (const auto range = field_range(content, 8); range.has_value()) {
        result.replace(range->first, range->second - range->first, formatted);
    } else {
        result += " 1 UNL ";
        result += formatted;
    }
    result += ending;
    return result;
}

[[nodiscard]] auto line_content(const std::string_view line) -> std::string_view {
    if (line.ends_with("\r\n")) {
        return line.substr(0, line.size() - 2);
    }
    if (line.ends_with('\n')) {
        return line.substr(0, line.size() - 1);
    }
    return line;
}

[[nodiscard]] auto line_ending(const std::string_view line) -> std::string_view {
    if (line.ends_with("\r\n")) {
        return "\r\n";
    }
    return line.ends_with('\n') ? "\n" : "";
}

[[nodiscard]] auto next_line(std::istream& input) -> std::optional<std::string> {
    std::string line;
    if (!std::getline(input, line)) {
        return std::nullopt;
    }
    line += '\n';
    return line;
}

auto write_record(const std::string_view record, const charges::ChargeAssignment* assignment,
                  std::ostream& output) -> void {
    if (assignment == nullptr) {
        std::print(output, "{}", record);
        return;
    }

    bool in_atom_section = false;
    std::size_t atom_index = 0;
    std::size_t start = 0;
    while (start < record.size()) {
        const auto newline = record.find('\n', start);
        const auto end = newline == std::string_view::npos ? record.size() : newline + 1;
        const auto line = record.substr(start, end - start);
        const auto content = line_content(line);

        if (content == atom_marker) {
            in_atom_section = true;
        } else if (content == bond_marker) {
            if (atom_index != assignment->charges.size()) {
                throw std::runtime_error{"MOL2 atom count does not match charge assignment"};
            }
            in_atom_section = false;
        }

        if (in_atom_section && !content.empty() && content.front() != '@') {
            if (atom_index == assignment->charges.size()) {
                throw std::runtime_error{"MOL2 atom count exceeds charge assignment"};
            }
            std::print(
                output, "{}",
                patch_atom_line(content, assignment->charges[atom_index++], line_ending(line)));
        } else {
            std::print(output, "{}", line);
        }
        start = end;
    }

    if (atom_index != assignment->charges.size()) {
        throw std::runtime_error{"MOL2 atom count does not match charge assignment"};
    }
}

auto write_preserving_records(std::istream& input, std::ostream& output,
                              const std::span<const charges::ChargeAssignment> assignments)
    -> void {
    std::size_t record_index = 0;
    auto pending = next_line(input);

    while (pending.has_value()) {
        if (line_content(*pending) != molecule_marker) {
            std::print(output, "{}", *pending);
            pending = next_line(input);
            continue;
        }

        auto record = std::move(*pending);
        pending.reset();
        while (const auto line = next_line(input)) {
            if (line_content(*line) == molecule_marker) {
                pending = *line;
                break;
            }
            record += *line;
        }

        if (record_index >= assignments.size()) {
            throw std::invalid_argument{"MOL2 source has more records than charge assignments"};
        }

        const auto& assignment = assignments[record_index];
        if (assignment.target.molecule_index != record_index) {
            throw std::invalid_argument{"MOL2 assignment order does not match source record order"};
        }
        write_record(record, std::addressof(assignment), output);
        ++record_index;
    }

    if (record_index != assignments.size()) {
        throw std::invalid_argument{"MOL2 source has fewer records than charge assignments"};
    }
}

} // namespace

Mol2Writer::Mol2Writer(std::ostream& output) : output_{std::addressof(output)} {}

auto Mol2Writer::write_preserving_source(
    const std::string& source_path,
    const std::span<const charges::ChargeAssignment> assignments) const -> void {
    auto input = common_output::open_source_file(source_path, "MOL2");
    write_preserving_records(input, *output_, assignments);
}

auto Mol2Writer::write_preserving_buffer(
    const std::string_view source,
    const std::span<const charges::ChargeAssignment> assignments) const -> void {
    auto input = std::istringstream{std::string{source}};
    write_preserving_records(input, *output_, assignments);
}

auto Mol2Writer::write_generated(const core::Molecule& molecule,
                                 const charges::ChargeAssignment& assignment) const -> void {
    const auto& conformer = common_output::assignment_conformer(molecule, assignment, "MOL2");

    std::print(*output_,
               "@<TRIPOS>MOLECULE\n{}\n{} {} 0 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n",
               molecule.name().empty() ? "chargefw" : molecule.name(), molecule.atom_count(),
               molecule.bond_count());
    for (std::size_t index = 0; index < molecule.atom_count(); ++index) {
        const auto& atom = molecule.atom(index);
        const auto& position = conformer[index];
        std::print(*output_, "{} {} {} {} {} {} 1 UNL {}\n", index + 1,
                   common_output::generated_atom_name(atom, index), position.x, position.y,
                   position.z, common_output::atom_element_symbol(atom),
                   common_output::formatted_charge(assignment.charges[index]));
    }

    std::print(*output_, "@<TRIPOS>BOND\n");
    for (std::size_t index = 0; index < molecule.bond_count(); ++index) {
        const auto& bond = molecule.bond(index);
        std::print(*output_, "{} {} {} {}\n", index + 1, bond.first_atom_index() + 1,
                   bond.second_atom_index() + 1, common_output::bond_type(bond.order()));
    }
}

} // namespace chargefw::adapters::native::mol2_output
