#include <chargefw/adapters/gemmi/mmcif_output.h>

#include <chargefw/core/periodic_table.h>

#include <gemmi/cif.hpp>
#include <gemmi/to_cif.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <format>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

namespace chargefw::adapters::gemmi::mmcif_output {
namespace {

constexpr auto metadata_category = "_sb_ncbr_partial_atomic_charges_meta.";
constexpr auto charges_category = "_sb_ncbr_partial_atomic_charges.";
constexpr auto dictionary_name = "mmcif_charges_v11.dic";
constexpr auto dictionary_version = "1.1";
constexpr auto dictionary_location =
    "https://sb-ncbr.github.io/charges-schema/schemas/mmcif_charges_v11.dic";

struct BlockMapping {
    std::vector<std::vector<std::string>> atom_site_ids;
    std::vector<std::string> model_ids;
};

auto ensure_dictionary(::gemmi::cif::Block& block) -> void;

[[nodiscard]] auto quote(const std::string_view value) -> std::string {
    return ::gemmi::cif::quote(std::string{value});
}

[[nodiscard]] auto block_name(const ImportedMoleculeRecord& record, const std::size_t index)
    -> std::string {
    auto candidate = record.identity.record_id.empty() ? std::string{record.molecule.name()}
                                                       : record.identity.record_id.display_string();
    if (candidate.empty()) {
        candidate = "molecule_" + std::to_string(index + 1);
    }
    for (auto& character : candidate) {
        if (std::isalnum(static_cast<unsigned char>(character)) == 0 && character != '_' &&
            character != '-' && character != '.') {
            character = '_';
        }
    }
    return candidate;
}

[[nodiscard]] auto unique_block_name(const ::gemmi::cif::Document& document, std::string candidate)
    -> std::string {
    if (document.find_block(candidate) == nullptr) {
        return candidate;
    }
    const auto base = candidate;
    for (std::size_t suffix = 2;; ++suffix) {
        candidate = base + "_" + std::to_string(suffix);
        if (document.find_block(candidate) == nullptr) {
            return candidate;
        }
    }
}

[[nodiscard]] auto atom_ids(const core::Molecule& molecule) -> std::vector<std::string> {
    std::vector<std::string> result;
    result.reserve(molecule.atom_count());
    std::unordered_set<std::string> used;
    for (std::size_t index = 0; index < molecule.atom_count(); ++index) {
        auto candidate = std::string{molecule.atom(index).name()};
        const auto valid =
            !candidate.empty() &&
            std::ranges::all_of(candidate, [](const unsigned char character) -> bool {
                return std::isalnum(character) != 0 || character == '_' || character == '-';
            });
        if (!valid || used.contains(candidate)) {
            candidate = std::string{core::element_symbol(molecule.atom(index).atomic_number())} +
                        std::to_string(index + 1);
        }
        while (used.contains(candidate)) {
            candidate += "_" + std::to_string(index + 1);
        }
        used.insert(candidate);
        result.push_back(std::move(candidate));
    }
    return result;
}

[[nodiscard]] auto bond_order(const core::BondOrder order) -> std::string {
    switch (order) {
    case core::BondOrder::SINGLE:
        return "sing";
    case core::BondOrder::DOUBLE:
        return "doub";
    case core::BondOrder::TRIPLE:
        return "trip";
    }
    throw std::runtime_error{"cannot write unsupported bond order to generated mmCIF"};
}

[[nodiscard]] auto round_trip_number(const double value) -> std::string {
    auto buffer = std::array<char, 64>{};
    const auto [end, error] =
        std::to_chars(buffer.data(), buffer.data() + buffer.size(), value,
                      std::chars_format::general, std::numeric_limits<double>::max_digits10);
    if (error != std::errc{}) {
        throw std::runtime_error{"cannot serialize floating-point mmCIF value"};
    }
    return std::string{buffer.data(), end};
}

[[nodiscard]] auto known_label(const std::optional<std::string>& value)
    -> std::optional<std::string_view> {
    if (!value.has_value() || value->empty() || *value == "." || *value == "?") {
        return std::nullopt;
    }
    return *value;
}

[[nodiscard]] auto label_or(const std::optional<std::string>& preferred,
                            const std::optional<std::string>& alternate,
                            const std::string_view fallback) -> std::string {
    if (const auto value = known_label(preferred); value.has_value()) {
        return std::string{*value};
    }
    if (const auto value = known_label(alternate); value.has_value()) {
        return std::string{*value};
    }
    return std::string{fallback};
}

struct OutputAssignments {
    const charges::ChargeAssignment* molecule = nullptr;
    std::vector<const charges::ChargeAssignment*> conformers;
};

[[nodiscard]] auto assignments_for_result(const ChargeCalculationResult& result,
                                          const std::size_t molecule_index) -> OutputAssignments {
    const auto& molecule = result.inputs()[molecule_index].molecule;
    auto assignments = OutputAssignments{
        .molecule = nullptr,
        .conformers = std::vector<const charges::ChargeAssignment*>(molecule.conformer_count())};
    for (const auto& assignment : result.execution().charges->assignments()) {
        if (assignment.target.molecule_index != molecule_index) {
            continue;
        }
        if (assignment.target.conformer_index.has_value()) {
            assignments.conformers[*assignment.target.conformer_index] = std::addressof(assignment);
        } else {
            assignments.molecule = std::addressof(assignment);
        }
    }
    return assignments;
}

[[nodiscard]] auto canonical_charge_site_id(const std::string_view value) -> bool {
    int id = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), id);
    return error == std::errc{} && end == value.data() + value.size() && id > 0 &&
           std::to_string(id) == value;
}

[[nodiscard]] auto attached_mapping(::gemmi::cif::Block& block,
                                    const ImportedMoleculeRecord& record) -> BlockMapping {
    if (!record.import_metadata.has_value() ||
        record.import_metadata->format != MolecularSourceFormat::mmcif) {
        throw std::invalid_argument{
            "charge attachment requires unchanged mmCIF-imported calculation inputs"};
    }
    if (record.identity.record_id.empty() ||
        record.identity.record_id.display_string() != block.name) {
        throw std::invalid_argument{
            "Gemmi target block identity does not match the calculation input"};
    }
    auto atom_sites = block.find("_atom_site.", {"id", "?pdbx_PDB_model_num"});
    if (atom_sites.length() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument{"Gemmi target contains too many atom sites"};
    }
    const auto& metadata = *record.import_metadata;
    auto result = BlockMapping{};
    result.atom_site_ids.reserve(metadata.conformers.size());
    result.model_ids.reserve(metadata.conformers.size());
    auto mapped_ids = std::unordered_set<std::string>{};
    for (const auto& conformer : metadata.conformers) {
        const auto model_id = conformer.id.value_or("1");
        result.model_ids.push_back(model_id);
        auto& ids = result.atom_site_ids.emplace_back();
        ids.reserve(conformer.sites.size());
        for (const auto& site : conformer.sites) {
            if (!site.id.has_value() || site.position >= atom_sites.length()) {
                throw std::invalid_argument{
                    "Gemmi target site mapping does not match the calculation input"};
            }
            const auto row_index = static_cast<int>(site.position);
            auto target_id = ::gemmi::cif::as_string(atom_sites[row_index][0]);
            if (target_id != *site.id || !mapped_ids.insert(target_id).second) {
                throw std::invalid_argument{
                    "Gemmi target site mapping does not match the calculation input"};
            }
            if (!canonical_charge_site_id(target_id)) {
                throw std::invalid_argument{
                    "Gemmi target atom IDs are not representable by the charge dictionary"};
            }
            if (atom_sites[row_index].has(1) &&
                ::gemmi::cif::as_string(atom_sites[row_index][1]) != model_id) {
                throw std::invalid_argument{
                    "Gemmi target model mapping does not match the calculation input"};
            }
            ids.push_back(std::move(target_id));
        }
    }
    return result;
}

[[nodiscard]] auto structural_labels(const ImportedMoleculeRecord& record,
                                     const std::size_t conformer_index,
                                     const std::size_t atom_index)
    -> const SourceStructuralLabels* {
    if (!record.import_metadata.has_value()) {
        return nullptr;
    }
    const auto& reference = record.import_metadata->conformers[conformer_index].sites[atom_index];
    return reference.structural_labels.has_value() ? std::addressof(*reference.structural_labels)
                                                   : nullptr;
}

auto validate_result_output(const ChargeCalculationResult& result) -> void {
    if (!result.execution().calculated() || !result.execution().charges.has_value()) {
        throw std::invalid_argument{"mmCIF output requires a successful calculation"};
    }
    if (result.inputs().empty()) {
        throw std::invalid_argument{"mmCIF output requires at least one molecule record"};
    }
    for (const auto& record : result.inputs()) {
        const auto& molecule = record.molecule;
        if (molecule.atom_count() == 0) {
            throw std::invalid_argument{"mmCIF output requires at least one atom per record"};
        }
        if (molecule.conformer_count() == 0) {
            throw std::invalid_argument{"mmCIF output requires coordinates"};
        }
        if (molecule.atom_count() != 0 &&
            molecule.conformer_count() >
                static_cast<std::size_t>(std::numeric_limits<int>::max()) / molecule.atom_count()) {
            throw std::invalid_argument{"mmCIF output contains too many atom sites"};
        }
        for (const auto& conformer : molecule.conformers()) {
            for (const auto& position : conformer.positions()) {
                if (!std::isfinite(position.x) || !std::isfinite(position.y) ||
                    !std::isfinite(position.z)) {
                    throw std::invalid_argument{"mmCIF output coordinates must be finite"};
                }
            }
        }
    }
    for (const auto& assignment : result.execution().charges->assignments()) {
        for (const auto charge : assignment.charges.values()) {
            if (charge < -5.0 || charge > 5.0) {
                throw std::invalid_argument{"mmCIF charge is outside the dictionary range [-5, 5]"};
            }
        }
    }
}

auto append_charge_metadata(::gemmi::cif::Table& metadata, const std::string_view id,
                            const charges::ChargeSet& charge_set,
                            const std::string_view generator_name,
                            const std::string_view generator_version) -> void {
    const auto parameter_set = charge_set.parameter_set_id();
    metadata.append_row({std::string{id}, "empirical", quote(charge_set.method_id()),
                         parameter_set.has_value() ? quote(*parameter_set) : ".",
                         quote(generator_name.empty() ? "unknown" : generator_name),
                         quote(generator_version.empty() ? "unknown" : generator_version)});
}

auto write_result_block(::gemmi::cif::Block& block, const ImportedMoleculeRecord& record,
                        const OutputAssignments& assignments, const charges::ChargeSet& charge_set,
                        const std::string_view generator_name,
                        const std::string_view generator_version) -> void {
    const auto& molecule = record.molecule;
    const auto generated_atom_ids = atom_ids(molecule);
    block.set_pair("_entry.id", quote(block.name));
    ensure_dictionary(block);

    block.init_loop("_atom_site.", {"id",
                                    "type_symbol",
                                    "label_atom_id",
                                    "label_alt_id",
                                    "label_comp_id",
                                    "label_asym_id",
                                    "label_entity_id",
                                    "label_seq_id",
                                    "pdbx_PDB_ins_code",
                                    "Cartn_x",
                                    "Cartn_y",
                                    "Cartn_z",
                                    "occupancy",
                                    "B_iso_or_equiv",
                                    "pdbx_formal_charge",
                                    "auth_seq_id",
                                    "auth_comp_id",
                                    "auth_asym_id",
                                    "auth_atom_id",
                                    "pdbx_PDB_model_num"});
    block.init_loop(metadata_category,
                    {"id", "type", "method", "parameter_set", "software_name", "software_version"});
    block.init_loop(charges_category, {"type_id", "atom_id", "charge"});

    auto atom_sites = block.find("_atom_site.", {"id",
                                                 "type_symbol",
                                                 "label_atom_id",
                                                 "label_alt_id",
                                                 "label_comp_id",
                                                 "label_asym_id",
                                                 "label_entity_id",
                                                 "label_seq_id",
                                                 "pdbx_PDB_ins_code",
                                                 "Cartn_x",
                                                 "Cartn_y",
                                                 "Cartn_z",
                                                 "occupancy",
                                                 "B_iso_or_equiv",
                                                 "pdbx_formal_charge",
                                                 "auth_seq_id",
                                                 "auth_comp_id",
                                                 "auth_asym_id",
                                                 "auth_atom_id",
                                                 "pdbx_PDB_model_num"});
    auto metadata = block.find(metadata_category, {"id", "type", "method", "parameter_set",
                                                   "software_name", "software_version"});
    auto charge_rows = block.find(charges_category, {"type_id", "atom_id", "charge"});

    auto molecule_type_id = std::optional<std::string>{};
    auto conformer_type_ids = std::vector<std::optional<std::string>>(molecule.conformer_count());
    auto next_type_id = std::size_t{1};
    if (assignments.molecule != nullptr) {
        molecule_type_id = std::to_string(next_type_id++);
        append_charge_metadata(metadata, *molecule_type_id, charge_set, generator_name,
                               generator_version);
    } else {
        for (std::size_t index = 0; index < assignments.conformers.size(); ++index) {
            conformer_type_ids[index] = std::to_string(next_type_id++);
            append_charge_metadata(metadata, *conformer_type_ids[index], charge_set, generator_name,
                                   generator_version);
        }
    }

    auto next_site_id = std::size_t{1};
    for (std::size_t conformer_index = 0; conformer_index < molecule.conformer_count();
         ++conformer_index) {
        const auto model_id = std::to_string(conformer_index + 1);
        const auto* assignment = assignments.molecule != nullptr
                                     ? assignments.molecule
                                     : assignments.conformers[conformer_index];
        const auto& type_id =
            molecule_type_id.has_value() ? *molecule_type_id : *conformer_type_ids[conformer_index];
        for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
            const auto& atom = molecule.atom(atom_index);
            const auto& position = molecule.conformer(conformer_index)[atom_index];
            const auto* labels = structural_labels(record, conformer_index, atom_index);
            const auto& fallback_atom = generated_atom_ids[atom_index];
            const auto label_atom =
                labels == nullptr
                    ? fallback_atom
                    : label_or(labels->label.atom, labels->author.atom, fallback_atom);
            const auto author_atom =
                labels == nullptr
                    ? fallback_atom
                    : label_or(labels->author.atom, labels->label.atom, fallback_atom);
            const auto label_residue =
                labels == nullptr ? std::string{"UNL"}
                                  : label_or(labels->label.residue, labels->author.residue, "UNL");
            const auto author_residue =
                labels == nullptr ? std::string{"UNL"}
                                  : label_or(labels->author.residue, labels->label.residue, "UNL");
            const auto label_chain = labels == nullptr
                                         ? std::string{"A"}
                                         : label_or(labels->label.chain, labels->author.chain, "A");
            const auto author_chain =
                labels == nullptr ? std::string{"A"}
                                  : label_or(labels->author.chain, labels->label.chain, "A");
            const auto label_sequence =
                labels == nullptr || !known_label(labels->label.sequence).has_value()
                    ? std::string{"."}
                    : std::string{*known_label(labels->label.sequence)};
            const auto author_sequence =
                labels == nullptr
                    ? std::string{"1"}
                    : known_label(labels->author.sequence)
                          .transform([](const auto value) { return std::string{value}; })
                          .value_or(".");
            const auto entity =
                labels == nullptr
                    ? std::string{"1"}
                    : known_label(labels->entity)
                          .transform([](const auto value) { return std::string{value}; })
                          .value_or("1");
            const auto alternate_location =
                labels == nullptr
                    ? std::string{"."}
                    : known_label(labels->alternate_location)
                          .transform([](const auto value) { return std::string{value}; })
                          .value_or(".");
            const auto insertion_code =
                labels == nullptr
                    ? std::string{"."}
                    : known_label(labels->insertion_code)
                          .transform([](const auto value) { return std::string{value}; })
                          .value_or(".");
            const auto site_id = std::to_string(next_site_id++);
            atom_sites.append_row({site_id,
                                   std::string{core::element_symbol(atom.atomic_number())},
                                   quote(label_atom),
                                   alternate_location == "." ? "." : quote(alternate_location),
                                   quote(label_residue),
                                   quote(label_chain),
                                   quote(entity),
                                   label_sequence == "." ? "." : quote(label_sequence),
                                   insertion_code == "." ? "." : quote(insertion_code),
                                   round_trip_number(position.x),
                                   round_trip_number(position.y),
                                   round_trip_number(position.z),
                                   "?",
                                   "?",
                                   std::to_string(atom.formal_charge()),
                                   author_sequence == "." ? "." : quote(author_sequence),
                                   quote(author_residue),
                                   quote(author_chain),
                                   quote(author_atom),
                                   model_id});
            charge_rows.append_row(
                {type_id, site_id, round_trip_number(assignment->charges[atom_index])});
        }
    }
}

[[nodiscard]] auto write_generated_block(::gemmi::cif::Block& block, const core::Molecule& molecule)
    -> BlockMapping {
    if (molecule.conformer_count() == 0) {
        throw std::runtime_error{"mmCIF output requires coordinates"};
    }

    const auto ids = atom_ids(molecule);
    block.set_pair("_entry.id", quote(block.name));

    auto& entity = block.init_loop("_entity.", {"id", "type", "src_method", "pdbx_description"});
    entity.add_row({"1", "non-polymer", "syn", quote(molecule.name())});
    auto& asym = block.init_loop("_struct_asym.", {"id", "entity_id"});
    asym.add_row({"A", "1"});
    auto& nonpoly = block.init_loop("_pdbx_entity_nonpoly.", {"entity_id", "name", "comp_id"});
    nonpoly.add_row({"1", quote(molecule.name()), "UNL"});

    auto& component =
        block.init_loop("_chem_comp.", {"id", "name", "type", "formula", "formula_weight"});
    component.add_row({"UNL", quote(molecule.name()), "NON-POLYMER", "?", "?"});

    auto& atom_types = block.init_loop("_atom_type.", {"symbol"});
    std::vector<std::string> atom_type_symbols;
    atom_type_symbols.reserve(molecule.atom_count());
    for (std::size_t index = 0; index < molecule.atom_count(); ++index) {
        auto symbol = std::string{core::element_symbol(molecule.atom(index).atomic_number())};
        if (std::ranges::find(atom_type_symbols, symbol) == atom_type_symbols.end()) {
            atom_type_symbols.push_back(std::move(symbol));
        }
    }
    for (const auto& symbol : atom_type_symbols) {
        atom_types.add_row({symbol});
    }

    auto& component_atoms =
        block.init_loop("_chem_comp_atom.", {"comp_id", "atom_id", "type_symbol", "charge"});
    for (std::size_t index = 0; index < molecule.atom_count(); ++index) {
        const auto& atom = molecule.atom(index);
        component_atoms.add_row({"UNL", quote(ids[index]),
                                 std::string{core::element_symbol(atom.atomic_number())},
                                 std::to_string(atom.formal_charge())});
    }

    if (molecule.bond_count() != 0) {
        auto& component_bonds = block.init_loop(
            "_chem_comp_bond.", {"comp_id", "atom_id_1", "atom_id_2", "value_order"});
        for (const auto& bond : molecule.bonds()) {
            component_bonds.add_row({"UNL", quote(ids[bond.first_atom_index()]),
                                     quote(ids[bond.second_atom_index()]),
                                     bond_order(bond.order())});
        }
    }

    auto& atom_sites = block.init_loop("_atom_site.", {"group_PDB",
                                                       "id",
                                                       "type_symbol",
                                                       "label_atom_id",
                                                       "label_alt_id",
                                                       "label_comp_id",
                                                       "label_asym_id",
                                                       "label_entity_id",
                                                       "label_seq_id",
                                                       "Cartn_x",
                                                       "Cartn_y",
                                                       "Cartn_z",
                                                       "occupancy",
                                                       "B_iso_or_equiv",
                                                       "pdbx_formal_charge",
                                                       "auth_seq_id",
                                                       "auth_comp_id",
                                                       "auth_asym_id",
                                                       "auth_atom_id",
                                                       "pdbx_PDB_model_num"});

    BlockMapping mapping;
    mapping.atom_site_ids.resize(molecule.conformer_count());
    mapping.model_ids.reserve(molecule.conformer_count());
    std::size_t site_id = 1;
    for (std::size_t conformer_index = 0; conformer_index < molecule.conformer_count();
         ++conformer_index) {
        const auto model_id = std::to_string(conformer_index + 1);
        mapping.model_ids.push_back(model_id);
        auto& site_ids = mapping.atom_site_ids[conformer_index];
        site_ids.reserve(molecule.atom_count());
        for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
            const auto& atom = molecule.atom(atom_index);
            const auto& position = molecule.conformer(conformer_index)[atom_index];
            const auto current_site_id = std::to_string(site_id++);
            site_ids.push_back(current_site_id);
            atom_sites.add_row({"HETATM",
                                current_site_id,
                                std::string{core::element_symbol(atom.atomic_number())},
                                quote(ids[atom_index]),
                                ".",
                                "UNL",
                                "A",
                                "1",
                                ".",
                                std::format("{:.6f}", position.x),
                                std::format("{:.6f}", position.y),
                                std::format("{:.6f}", position.z),
                                "1.0",
                                "0.0",
                                std::to_string(atom.formal_charge()),
                                "1",
                                "UNL",
                                "A",
                                quote(ids[atom_index]),
                                quote(model_id)});
        }
    }
    return mapping;
}

auto erase_category(::gemmi::cif::Block& block, const std::string_view category) -> void {
    if (block.has_mmcif_category(std::string{category})) {
        block.find_mmcif_category(std::string{category}).erase();
    }
}

auto ensure_dictionary(::gemmi::cif::Block& block) -> void {
    auto table =
        block.find_or_add("_audit_conform.", {"dict_name", "dict_version", "dict_location"});
    for (auto row : table) {
        if (::gemmi::cif::as_string(row[0]) == dictionary_name) {
            return;
        }
    }
    table.ensure_loop();
    table.append_row({dictionary_name, dictionary_version, dictionary_location});
}

auto write_charges(::gemmi::cif::Block& block, const BlockMapping& mapping,
                   const core::Molecule& molecule,
                   const std::span<const charges::ChargeAssignment> assignments,
                   const charges::ChargeSet& charge_set, const std::string_view generator_name,
                   const std::string_view generator_version) -> void {
    erase_category(block, metadata_category);
    erase_category(block, charges_category);
    ensure_dictionary(block);
    auto assignment_id = std::size_t{1};
    block
        .find_or_add(metadata_category,
                     {"id", "type", "method", "parameter_set", "software_name", "software_version"})
        .ensure_loop();
    block.find_or_add(charges_category, {"type_id", "atom_id", "charge"}).ensure_loop();
    auto metadata = block.find(metadata_category, {"id", "type", "method", "parameter_set",
                                                   "software_name", "software_version"});
    auto charge_rows = block.find(charges_category, {"type_id", "atom_id", "charge"});

    for (const auto& assignment : assignments) {
        if (assignment.charges.size() != molecule.atom_count()) {
            throw std::runtime_error{"charge assignment size does not match mmCIF molecule"};
        }
        if (assignment.target.conformer_index.has_value() &&
            *assignment.target.conformer_index >= mapping.atom_site_ids.size()) {
            throw std::runtime_error{"charge assignment conformer is missing from mmCIF block"};
        }
        const auto id = std::to_string(assignment_id++);
        const auto parameter_set = charge_set.parameter_set_id();
        metadata.append_row({id, "empirical", quote(charge_set.method_id()),
                             parameter_set.has_value() ? quote(*parameter_set) : ".",
                             quote(generator_name.empty() ? "unknown" : generator_name),
                             quote(generator_version.empty() ? "unknown" : generator_version)});
        const auto first_mapping = assignment.target.conformer_index.value_or(0);
        const auto mapping_count = assignment.target.conformer_index.has_value()
                                       ? first_mapping + 1
                                       : mapping.atom_site_ids.size();
        for (std::size_t mapping_index = first_mapping; mapping_index < mapping_count;
             ++mapping_index) {
            for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
                charge_rows.append_row({id, quote(mapping.atom_site_ids[mapping_index][atom_index]),
                                        round_trip_number(assignment.charges[atom_index])});
            }
        }
    }
}

[[nodiscard]] auto assignments_for(const charges::ChargeSet& charge_set,
                                   const std::size_t molecule_index)
    -> std::vector<charges::ChargeAssignment> {
    std::vector<charges::ChargeAssignment> result;
    for (const auto& assignment : charge_set.assignments()) {
        if (assignment.target.molecule_index == molecule_index) {
            result.push_back(assignment);
        }
    }
    if (result.empty()) {
        throw std::runtime_error{"no charge assignments for mmCIF output molecule " +
                                 std::to_string(molecule_index + 1)};
    }
    return result;
}

} // namespace

MmcifWriter::MmcifWriter(std::ostream& output) : output_{std::addressof(output)} {}

auto MmcifWriter::write(const ChargeCalculationResult& result,
                        const std::string_view generator_name,
                        const std::string_view generator_version) const -> void {
    validate_result_output(result);
    auto document = ::gemmi::cif::Document{};
    for (std::size_t index = 0; index < result.inputs().size(); ++index) {
        const auto& record = result.inputs()[index];
        auto& block =
            document.add_new_block(unique_block_name(document, block_name(record, index)));
        write_result_block(block, record, assignments_for_result(result, index),
                           *result.execution().charges, generator_name, generator_version);
    }
    ::gemmi::cif::write_cif_to_stream(*output_, document);
    if (!*output_) {
        throw std::runtime_error{"failed to write mmCIF output"};
    }
}

auto MmcifWriter::write_attached(const ChargeCalculationResult& result,
                                 const ::gemmi::cif::Document& source, const bool overwrite,
                                 const std::string_view generator_name,
                                 const std::string_view generator_version) const -> void {
    validate_result_output(result);
    auto document = source;
    auto blocks = std::vector<::gemmi::cif::Block*>{};
    for (auto& block : document.blocks) {
        if (block.has_mmcif_category("_atom_site.")) {
            blocks.push_back(std::addressof(block));
        }
    }
    if (blocks.size() != result.inputs().size()) {
        throw std::invalid_argument{
            "Gemmi target molecule count does not match the calculation input"};
    }
    for (std::size_t index = 0; index < blocks.size(); ++index) {
        auto& block = *blocks[index];
        if (!overwrite && (block.has_mmcif_category(metadata_category) ||
                           block.has_mmcif_category(charges_category))) {
            throw std::invalid_argument{"Gemmi target already contains partial charge categories"};
        }
        const auto mapping = attached_mapping(block, result.inputs()[index]);
        write_charges(block, mapping, result.inputs()[index].molecule,
                      assignments_for(*result.execution().charges, index),
                      *result.execution().charges, generator_name, generator_version);
    }
    ::gemmi::cif::write_cif_to_stream(*output_, document);
    if (!*output_) {
        throw std::runtime_error{"failed to write mmCIF output"};
    }
}

auto MmcifWriter::write_generated(const std::span<const ImportedMoleculeRecord> records,
                                  const charges::ChargeSet& charge_set,
                                  const std::string_view generator_name,
                                  const std::string_view generator_version) const -> void {
    if (records.empty()) {
        throw std::invalid_argument{"mmCIF output requires at least one molecule record"};
    }

    ::gemmi::cif::Document document;
    for (std::size_t record_index = 0; record_index < records.size(); ++record_index) {
        const auto& record = records[record_index];
        const auto assignments = assignments_for(charge_set, record_index);
        auto& block =
            document.add_new_block(unique_block_name(document, block_name(record, record_index)));
        const auto mapping = write_generated_block(block, record.molecule);
        write_charges(block, mapping, record.molecule, assignments, charge_set, generator_name,
                      generator_version);
    }

    ::gemmi::cif::write_cif_to_stream(*output_, document);
    if (!*output_) {
        throw std::runtime_error{"failed to write mmCIF output"};
    }
}

} // namespace chargefw::adapters::gemmi::mmcif_output
