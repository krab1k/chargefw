#include <chargefw/adapters/native/mol2_output.h>

#include "common_output.h"

#include <cstddef>
#include <fstream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::adapters::native::mol2_output {
namespace {

constexpr std::string_view molecule_marker{"@<TRIPOS>MOLECULE"};
constexpr std::string_view atom_marker{"@<TRIPOS>ATOM"};
constexpr std::string_view bond_marker{"@<TRIPOS>BOND"};
[[nodiscard]] auto patch_atom_line(const std::string_view content, const double charge,
                                   const std::string_view ending) -> std::string {
    auto fields = std::vector<std::string>{};
    auto input = std::istringstream{std::string{content}};
    for (auto field = std::string{}; input >> field;) {
        fields.push_back(std::move(field));
    }
    if (fields.size() < 6) {
        throw std::runtime_error{"invalid MOL2 atom record while writing charges"};
    }

    if (fields.size() >= 9) {
        fields[8] = common_output::formatted_charge(charge);
    } else {
        while (fields.size() < 8) {
            fields.emplace_back(fields.size() == 6 ? "1" : "CHARGEFW");
        }
        fields.push_back(common_output::formatted_charge(charge));
    }

    auto output = std::ostringstream{};
    for (std::size_t index = 0; index < fields.size(); ++index) {
        if (index != 0) {
            output << ' ';
        }
        output << fields[index];
    }
    output << ending;
    return output.str();
}

} // namespace

Mol2Writer::Mol2Writer(std::ostream& output) : output_{std::addressof(output)} {}

auto Mol2Writer::write_preserving_source(const std::string& source_path,
                                         const std::size_t record_index,
                                         const charges::ChargeAssignment& assignment) const
    -> void {
    auto input = std::ifstream{source_path, std::ios::binary};
    if (!input) {
        throw std::runtime_error{"unable to open MOL2 source file: " + source_path};
    }

    std::size_t current_record = 0;
    bool selected_record = false;
    bool in_atom_section = false;
    std::size_t atom_index = 0;
    const auto source = std::string{std::istreambuf_iterator<char>{input}, {}};
    std::size_t line_start = 0;
    while (line_start < source.size()) {
        const auto line_end = source.find('\n', line_start);
        const auto has_newline = line_end != std::string::npos;
        const auto next_line_start = has_newline ? line_end + 1 : source.size();
        const auto line_size = next_line_start - line_start;
        const auto line = std::string_view{source}.substr(line_start, line_size);
        auto ending = std::string_view{};
        if (has_newline) {
            ending = line.ends_with("\r\n") ? "\r\n" : "\n";
        }
        const auto content = line.substr(0, line.size() - std::string_view{ending}.size());

        if (content == molecule_marker) {
            if (selected_record && atom_index != assignment.charges.size()) {
                throw std::runtime_error{"MOL2 atom count does not match charge assignment"};
            }
            selected_record = current_record++ == record_index;
            in_atom_section = false;
            atom_index = 0;
        } else if (selected_record && content == atom_marker) {
            in_atom_section = true;
        } else if (selected_record && content == bond_marker) {
            if (atom_index != assignment.charges.size()) {
                throw std::runtime_error{"MOL2 atom count does not match charge assignment"};
            }
            in_atom_section = false;
        }

        if (selected_record && in_atom_section && !content.empty() && content.front() != '@') {
            if (atom_index == assignment.charges.size()) {
                throw std::runtime_error{"MOL2 atom count exceeds charge assignment"};
            }
            *output_ << patch_atom_line(content, assignment.charges.at(atom_index++), ending);
        } else {
            *output_ << line;
        }
        line_start = next_line_start;
    }

    if (!selected_record) {
        throw std::runtime_error{"MOL2 record index not found in source file"};
    }
    if (atom_index != assignment.charges.size()) {
        throw std::runtime_error{"MOL2 atom count does not match charge assignment"};
    }
}

auto Mol2Writer::write_generated(const core::Molecule& molecule,
                                 const charges::ChargeAssignment& assignment) const -> void {
    common_output::validate_assignment(assignment.charges, molecule.atom_count());
    if (!assignment.target.conformer_index.has_value()) {
        throw std::invalid_argument{
            "generated MOL2 output requires a conformer-specific assignment"};
    }
    const auto conformer_index = *assignment.target.conformer_index;
    if (conformer_index >= molecule.conformer_count()) {
        throw std::invalid_argument{"charge assignment references an unavailable conformer"};
    }

    const auto& conformer = molecule.conformer(conformer_index);
    if (conformer.size() != molecule.atom_count()) {
        throw std::invalid_argument{
            "MOL2 conformer coordinate count does not match molecule atom count"};
    }

    *output_ << "@<TRIPOS>MOLECULE\n"
             << (molecule.name().empty() ? "chargefw" : molecule.name()) << '\n'
             << molecule.atom_count() << ' ' << molecule.bond_count() << " 0 0 0\nSMALL\n"
             << "USER_CHARGES\n\n@<TRIPOS>ATOM\n";
    for (std::size_t index = 0; index < molecule.atom_count(); ++index) {
        const auto& atom = molecule.atom(index);
        const auto& position = conformer[index];
        *output_ << index + 1 << ' ' << common_output::generated_atom_name(atom, index) << ' '
                 << position.x << ' ' << position.y << ' ' << position.z << ' '
                 << common_output::mol2_atom_type(atom) << " 1 CHARGEFW "
                 << common_output::formatted_charge(assignment.charges[index]) << '\n';
    }

    *output_ << "@<TRIPOS>BOND\n";
    for (std::size_t index = 0; index < molecule.bond_count(); ++index) {
        const auto& bond = molecule.bond(index);
        *output_ << index + 1 << ' ' << bond.first_atom_index() + 1 << ' '
                 << bond.second_atom_index() + 1 << ' '
                 << common_output::mol2_bond_type(bond.order()) << '\n';
    }
}

} // namespace chargefw::adapters::native::mol2_output
