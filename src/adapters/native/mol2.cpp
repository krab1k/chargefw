#include <chargefw/adapters/native/mol2.h>

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/core/position.h>

#include <charconv>
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

namespace chargefw::adapters::native::mol2 {
namespace {

constexpr std::string_view molecule_marker{"@<TRIPOS>MOLECULE"};
constexpr std::string_view atom_marker{"@<TRIPOS>ATOM"};
constexpr std::string_view bond_marker{"@<TRIPOS>BOND"};

[[nodiscard]] auto read_line(std::istream& input, std::size_t& line) -> std::string {
    std::string result;

    if (!std::getline(input, result)) {
        throw std::runtime_error{"unexpected end of MOL2 record"};
    }

    ++line;
    return result;
}

[[nodiscard]] auto parse_int(const std::string_view value, const std::string_view field) -> int {
    int result = 0;
    const auto [end, error] = std::from_chars(value.begin(), value.end(), result);

    if (error != std::errc{} || end != value.end()) {
        throw std::runtime_error{"invalid " + std::string{field}};
    }

    return result;
}

[[nodiscard]] auto parse_double(const std::string_view value, const std::string_view field)
    -> double {
    std::istringstream input{std::string{value}};
    double result = 0.0;
    char remaining = '\0';

    if (!(input >> result) || (input >> remaining)) {
        throw std::runtime_error{"invalid " + std::string{field}};
    }

    return result;
}

[[nodiscard]] auto element_symbol(const std::string_view atom_type) -> std::string_view {
    const auto separator = atom_type.find('.');
    const auto symbol = atom_type.substr(0, separator);

    if (symbol.empty() || symbol == "Du" || symbol == "LP" || symbol == "Any" || symbol == "Hal" ||
        symbol == "Het" || symbol == "Hev") {
        throw std::runtime_error{"unsupported MOL2 atom type '" + std::string{atom_type} + "'"};
    }

    return symbol;
}

[[nodiscard]] auto bond_order(const std::string_view type) -> core::BondOrder {
    if (type == "1") {
        return core::BondOrder::SINGLE;
    }

    if (type == "2") {
        return core::BondOrder::DOUBLE;
    }

    if (type == "3") {
        return core::BondOrder::TRIPLE;
    }

    if (type == "ar") {
        return core::BondOrder::AROMATIC;
    }

    throw std::runtime_error{"unsupported MOL2 bond type '" + std::string{type} + "'"};
}

auto consume_to_next_molecule(std::istream& input) -> std::optional<std::string> {
    std::string line;

    while (std::getline(input, line)) {
        if (line == molecule_marker) {
            return line;
        }
    }

    return std::nullopt;
}

[[nodiscard]] auto parse_record(std::istream& input, MoleculeRecordIdentity identity)
    -> ::chargefw::adapters::MoleculeRecordResult {
    std::size_t line = 1;

    try {
        const auto name = read_line(input, line);
        const auto counts_line = read_line(input, line);
        std::istringstream counts{counts_line};
        std::string atom_count_text;
        std::string bond_count_text;

        if (!(counts >> atom_count_text >> bond_count_text)) {
            throw std::runtime_error{"invalid MOL2 counts"};
        }

        const auto atom_count = parse_int(atom_count_text, "MOL2 atom count");
        const auto bond_count = parse_int(bond_count_text, "MOL2 bond count");

        if (atom_count <= 0 || bond_count < 0) {
            throw std::runtime_error{"invalid MOL2 counts"};
        }

        while (read_line(input, line) != atom_marker) {
        }

        std::vector<core::Atom> atoms;
        std::vector<core::Position> positions;
        std::unordered_map<int, std::size_t> atom_indices;
        bool has_partial_charges = false;
        atoms.reserve(static_cast<std::size_t>(atom_count));
        positions.reserve(static_cast<std::size_t>(atom_count));

        for (int index = 0; index < atom_count; ++index) {
            std::istringstream atom_line{read_line(input, line)};
            std::string source_id_text;
            std::string atom_name;
            std::string x_text;
            std::string y_text;
            std::string z_text;
            std::string atom_type;
            std::string subst_id;
            std::string subst_name;
            std::string charge_text;

            if (!(atom_line >> source_id_text >> atom_name >> x_text >> y_text >> z_text >>
                  atom_type)) {
                throw std::runtime_error{"invalid MOL2 atom record"};
            }

            static_cast<void>(atom_line >> subst_id >> subst_name >> charge_text);
            const auto source_id = parse_int(source_id_text, "MOL2 atom id");

            if (source_id <= 0 || atom_indices.contains(source_id)) {
                throw std::runtime_error{"duplicate or invalid MOL2 atom id"};
            }

            if (!charge_text.empty()) {
                const auto charge = parse_double(charge_text, "MOL2 atom charge");

                if (charge != 0.0) {
                    has_partial_charges = true;
                }
            }

            const auto& element = core::periodic_table().element(element_symbol(atom_type));
            atom_indices.emplace(source_id, atoms.size());
            atoms.emplace_back(element.atomic_number, 0, std::move(atom_name));
            positions.push_back(core::Position{.x = parse_double(x_text, "MOL2 x coordinate"),
                                               .y = parse_double(y_text, "MOL2 y coordinate"),
                                               .z = parse_double(z_text, "MOL2 z coordinate")});
        }

        if (read_line(input, line) != bond_marker) {
            throw std::runtime_error{"expected MOL2 BOND section"};
        }

        std::vector<core::Bond> bonds;
        bonds.reserve(static_cast<std::size_t>(bond_count));

        for (int index = 0; index < bond_count; ++index) {
            std::istringstream bond_line{read_line(input, line)};
            std::string bond_id;
            std::string first_text;
            std::string second_text;
            std::string type;

            if (!(bond_line >> bond_id >> first_text >> second_text >> type)) {
                throw std::runtime_error{"invalid MOL2 bond record"};
            }

            const auto first_source_id = parse_int(first_text, "MOL2 bond atom");
            const auto second_source_id = parse_int(second_text, "MOL2 bond atom");
            const auto first = atom_indices.find(first_source_id);
            const auto second = atom_indices.find(second_source_id);

            if (first == atom_indices.end() || second == atom_indices.end()) {
                throw std::runtime_error{"MOL2 bond references an unknown atom"};
            }

            bonds.emplace_back(first->second, second->second, bond_order(type));
        }

        if (identity.record_id.empty()) {
            identity.record_id = name;
        }

        auto atom_mapping = std::vector<std::optional<std::size_t>>{};
        atom_mapping.reserve(atoms.size());

        for (std::size_t index = 0; index < atoms.size(); ++index) {
            atom_mapping.emplace_back(index);
        }

        auto diagnostics = std::vector<MoleculeRecordDiagnostic>{};

        if (has_partial_charges) {
            diagnostics.push_back(MoleculeRecordDiagnostic{
                .message = "MOL2 partial charges were ignored; they are not formal charges and "
                           "were not used for calculation.",
                .line = std::nullopt});
        }

        return ImportedMoleculeRecord{
            .molecule = core::Molecule{std::move(atoms),
                                       std::move(bonds),
                                       {core::Conformer{std::move(positions), "input"}},
                                       identity.record_id},
            .identity = std::move(identity),
            .mapping = {.atom_indices = std::move(atom_mapping), .conformer_indices = {0}},
            .diagnostics = std::move(diagnostics),
            .source = std::nullopt};
    } catch (const std::exception& error) {
        return std::unexpected(MoleculeRecordError{
            .identity = std::move(identity), .message = error.what(), .line = std::optional{line}});
    }
}

} // namespace

Mol2Reader::Mol2Reader(std::istream& input, std::string source)
    : input_{std::addressof(input)}, source_{std::move(source)} {}

auto Mol2Reader::next() -> std::optional<::chargefw::adapters::MoleculeRecordResult> {
    if (!consume_to_next_molecule(*input_).has_value()) {
        return std::nullopt;
    }

    const auto identity =
        MoleculeRecordIdentity{.source = source_, .record_index = record_index_++, .record_id = {}};
    return parse_record(*input_, identity);
}

} // namespace chargefw::adapters::native::mol2
