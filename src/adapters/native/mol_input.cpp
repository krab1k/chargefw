#include <chargefw/adapters/native/mol_input.h>

#include "common_input.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/core/position.h>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <expected>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace chargefw::adapters::native::mol_input {
namespace {

namespace common = chargefw::adapters::native::common_input;

constexpr std::string_view v30_prefix{"M  V30 "};

[[nodiscard]] auto atom_from_symbol(std::string_view symbol, const int formal_charge)
    -> core::Atom {
    const auto normalized = common::trim(symbol);

    if (normalized.empty() || normalized == "*" || normalized == "A" || normalized == "Q" ||
        normalized == "L" || normalized == "LP" || normalized == "R" || normalized == "R#") {
        throw std::runtime_error{"unsupported query atom '" + std::string{normalized} + "'"};
    }

    const auto& element = core::periodic_table().element(normalized);
    return core::Atom{element.atomic_number, formal_charge, std::string{normalized}};
}

struct ParsedMolecule {
    std::vector<core::Atom> atoms;
    std::vector<core::Bond> bonds;
    std::vector<core::Position> positions;
    std::vector<MoleculeRecordDiagnostic> diagnostics;
};

auto parse_v2000(std::istream& input, const std::string_view counts, std::size_t& line)
    -> ParsedMolecule {
    const auto atom_count = common::parse_trimmed_int(
        common::fixed_field(counts, 0, 3, "V2000 atom count"), "V2000 atom count");
    const auto bond_count = common::parse_trimmed_int(
        common::fixed_field(counts, 3, 3, "V2000 bond count"), "V2000 bond count");

    if (atom_count <= 0 || bond_count < 0) {
        throw std::runtime_error{"invalid V2000 counts"};
    }

    auto result = ParsedMolecule{};
    result.atoms.reserve(static_cast<std::size_t>(atom_count));
    result.positions.reserve(static_cast<std::size_t>(atom_count));
    result.bonds.reserve(static_cast<std::size_t>(bond_count));

    for (int index = 0; index < atom_count; ++index) {
        const auto atom_line = common::read_line(input, line, "MOL");
        const auto x = common::parse_trimmed_double(
            common::fixed_field(atom_line, 0, 10, "V2000 x coordinate"), "V2000 x coordinate");
        const auto y = common::parse_trimmed_double(
            common::fixed_field(atom_line, 10, 10, "V2000 y coordinate"), "V2000 y coordinate");
        const auto z = common::parse_trimmed_double(
            common::fixed_field(atom_line, 20, 10, "V2000 z coordinate"), "V2000 z coordinate");
        const auto symbol = common::fixed_field(atom_line, 31, 3, "V2000 atom symbol");

        result.atoms.push_back(atom_from_symbol(symbol, 0));
        result.positions.push_back(core::Position{.x = x, .y = y, .z = z});
    }

    for (int index = 0; index < bond_count; ++index) {
        const auto bond_line = common::read_line(input, line, "MOL");
        const auto first = common::parse_trimmed_int(
            common::fixed_field(bond_line, 0, 3, "V2000 bond atom"), "V2000 bond atom");
        const auto second = common::parse_trimmed_int(
            common::fixed_field(bond_line, 3, 3, "V2000 bond atom"), "V2000 bond atom");
        const auto order = common::parse_trimmed_int(
            common::fixed_field(bond_line, 6, 3, "V2000 bond order"), "V2000 bond order");

        if (first <= 0 || second <= 0 || first > atom_count || second > atom_count) {
            throw std::runtime_error{"V2000 bond references an unknown atom"};
        }

        result.bonds.emplace_back(static_cast<std::size_t>(first - 1),
                                  static_cast<std::size_t>(second - 1),
                                  common::numeric_bond_order(order));
    }

    while (true) {
        const auto property_line = common::read_line(input, line, "MOL");

        if (property_line == "M  END") {
            return result;
        }

        if (!property_line.starts_with("M  ")) {
            throw std::runtime_error{"unexpected V2000 property line"};
        }

        if (!property_line.starts_with("M  CHG")) {
            throw std::runtime_error{"unsupported V2000 property '" + property_line.substr(0, 6) +
                                     "'"};
        }

        const auto count = common::parse_trimmed_int(
            common::fixed_field(property_line, 6, 3, "M  CHG count"), "M  CHG count");

        if (count < 0 || property_line.size() < 9 + static_cast<std::size_t>(count) * 8) {
            throw std::runtime_error{"invalid M  CHG record"};
        }

        for (int index = 0; index < count; ++index) {
            const auto offset = 9 + static_cast<std::size_t>(index) * 8;
            const auto atom_number = common::parse_trimmed_int(
                common::fixed_field(property_line, offset, 4, "M  CHG atom"), "M  CHG atom");
            const auto charge = common::parse_trimmed_int(
                common::fixed_field(property_line, offset + 4, 4, "M  CHG charge"),
                "M  CHG charge");

            if (atom_number <= 0 || atom_number > atom_count) {
                throw std::runtime_error{"M  CHG references an unknown atom"};
            }

            const auto atom_index = static_cast<std::size_t>(atom_number - 1);
            const auto& atom = result.atoms[atom_index];
            result.atoms[atom_index] =
                core::Atom{atom.atomic_number(), charge, std::string{atom.name()}};
        }
    }
}

[[nodiscard]] auto v30_payload(std::string_view line) -> std::string_view {
    if (!line.starts_with(v30_prefix)) {
        throw std::runtime_error{"expected V3000 record"};
    }

    return line.substr(v30_prefix.size());
}

[[nodiscard]] auto read_v30_line(std::istream& input, std::size_t& line) -> std::string {
    auto result = common::read_line(input, line, "MOL");

    if (!result.starts_with(v30_prefix)) {
        throw std::runtime_error{"expected V3000 record"};
    }

    while (!result.empty() && result.back() == '-') {
        result.pop_back();
        const auto continuation = common::read_line(input, line, "MOL");

        if (!continuation.starts_with(v30_prefix)) {
            throw std::runtime_error{"invalid V3000 continuation"};
        }

        result += continuation.substr(v30_prefix.size());
    }

    return result;
}

[[nodiscard]] auto v30_tokens(std::string_view line) -> std::vector<std::string> {
    std::istringstream input{std::string{v30_payload(line)}};
    std::vector<std::string> result;
    std::string token;

    while (input >> token) {
        result.push_back(std::move(token));
    }

    return result;
}

auto parse_v3000(std::istream& input, std::size_t& line) -> ParsedMolecule {
    if (v30_payload(read_v30_line(input, line)) != "BEGIN CTAB") {
        throw std::runtime_error{"expected V3000 BEGIN CTAB"};
    }

    const auto counts = v30_tokens(read_v30_line(input, line));

    if (counts.size() < 3 || counts[0] != "COUNTS") {
        throw std::runtime_error{"expected V3000 COUNTS"};
    }

    const auto atom_count = common::parse_int(counts[1], "V3000 atom count");
    const auto bond_count = common::parse_int(counts[2], "V3000 bond count");

    if (atom_count <= 0 || bond_count < 0) {
        throw std::runtime_error{"invalid V3000 counts"};
    }

    if (v30_payload(read_v30_line(input, line)) != "BEGIN ATOM") {
        throw std::runtime_error{"expected V3000 BEGIN ATOM"};
    }

    auto result = ParsedMolecule{};
    result.atoms.reserve(static_cast<std::size_t>(atom_count));
    result.positions.reserve(static_cast<std::size_t>(atom_count));
    result.bonds.reserve(static_cast<std::size_t>(bond_count));
    std::unordered_map<int, std::size_t> atom_indices;

    for (int index = 0; index < atom_count; ++index) {
        const auto atom = v30_tokens(read_v30_line(input, line));

        if (atom.size() < 6) {
            throw std::runtime_error{"invalid V3000 atom record"};
        }

        const auto source_id = common::parse_int(atom[0], "V3000 atom id");

        if (source_id <= 0 || atom_indices.contains(source_id)) {
            throw std::runtime_error{"duplicate or invalid V3000 atom id"};
        }

        int formal_charge = 0;

        for (std::size_t attribute_index = 6; attribute_index < atom.size(); ++attribute_index) {
            const auto& attribute = atom[attribute_index];

            if (attribute.starts_with("CHG=")) {
                formal_charge =
                    common::parse_int(std::string_view{attribute}.substr(4), "V3000 CHG value");
            } else if (attribute.starts_with("CFG=")) {
                result.diagnostics.push_back(MoleculeRecordDiagnostic{
                    .message = "V3000 CFG stereochemical attribute was ignored", .line = line});
            } else {
                throw std::runtime_error{"unsupported V3000 atom attribute '" + attribute + "'"};
            }
        }

        atom_indices.emplace(source_id, result.atoms.size());
        result.atoms.push_back(atom_from_symbol(atom[1], formal_charge));
        result.positions.push_back(
            core::Position{.x = common::parse_double(atom[2], "V3000 x coordinate"),
                           .y = common::parse_double(atom[3], "V3000 y coordinate"),
                           .z = common::parse_double(atom[4], "V3000 z coordinate")});
    }

    if (v30_payload(read_v30_line(input, line)) != "END ATOM" ||
        v30_payload(read_v30_line(input, line)) != "BEGIN BOND") {
        throw std::runtime_error{"expected V3000 bond block"};
    }

    for (int index = 0; index < bond_count; ++index) {
        const auto bond = v30_tokens(read_v30_line(input, line));

        if (bond.size() != 4) {
            throw std::runtime_error{"unsupported V3000 bond attributes"};
        }

        const auto first_source_id = common::parse_int(bond[2], "V3000 bond atom");
        const auto second_source_id = common::parse_int(bond[3], "V3000 bond atom");
        const auto first = atom_indices.find(first_source_id);
        const auto second = atom_indices.find(second_source_id);

        if (first == atom_indices.end() || second == atom_indices.end()) {
            throw std::runtime_error{"V3000 bond references an unknown atom"};
        }

        result.bonds.emplace_back(
            first->second, second->second,
            common::numeric_bond_order(common::parse_int(bond[1], "V3000 bond order")));
    }

    if (v30_payload(read_v30_line(input, line)) != "END BOND" ||
        v30_payload(read_v30_line(input, line)) != "END CTAB" ||
        common::read_line(input, line, "MOL") != "M  END") {
        throw std::runtime_error{"invalid V3000 CTAB termination"};
    }

    return result;
}

[[nodiscard]] auto make_record(ParsedMolecule parsed, MoleculeRecordIdentity identity)
    -> ImportedMoleculeRecord {
    return common::make_record(std::move(parsed.atoms), std::move(parsed.bonds),
                               {core::Conformer{std::move(parsed.positions), "input"}},
                               std::move(identity), {}, std::move(parsed.diagnostics));
}

} // namespace

auto parse_mol(std::istream& input, MoleculeRecordIdentity identity)
    -> ::chargefw::adapters::MoleculeRecordResult {
    std::size_t line = 0;

    try {
        const auto name = common::read_line(input, line, "MOL");
        static_cast<void>(common::read_line(input, line, "MOL"));
        static_cast<void>(common::read_line(input, line, "MOL"));
        const auto counts = common::read_line(input, line, "MOL");

        if (identity.record_id.empty()) {
            identity.record_id = name;
        }

        ParsedMolecule parsed;

        if (counts.contains("V2000")) {
            parsed = parse_v2000(input, counts, line);
        } else if (counts.contains("V3000")) {
            parsed = parse_v3000(input, line);
        } else {
            throw std::runtime_error{"unsupported MOL version"};
        }

        return make_record(std::move(parsed), std::move(identity));
    } catch (const std::exception& error) {
        return std::unexpected(
            MoleculeRecordError{.identity = std::move(identity),
                                .message = error.what(),
                                .line = line == 0 ? std::nullopt : std::optional{line}});
    }
}

MolReader::MolReader(std::istream& input, std::string source)
    : input_{std::addressof(input)}, source_{std::move(source)} {}

auto MolReader::next() -> std::optional<::chargefw::adapters::MoleculeRecordResult> {
    if (consumed_ || input_->peek() == std::char_traits<char>::eof()) {
        return std::nullopt;
    }

    consumed_ = true;
    return parse_mol(*input_, {.source = source_, .record_index = 0, .record_id = {}});
}

} // namespace chargefw::adapters::native::mol_input
