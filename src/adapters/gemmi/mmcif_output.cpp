#include <chargefw/adapters/gemmi/mmcif_output.h>

#include "selection.h"

#include <chargefw/core/periodic_table.h>

#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <format>
#include <memory>
#include <optional>
#include <ostream>
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
constexpr auto dictionary_name = "mmcif_charges_v10.dic";
constexpr auto dictionary_version = "1.0";
constexpr auto dictionary_location =
    "https://sb-ncbr.github.io/charges-schema/schemas/mmcif_charges_v10.dic";

struct BlockMapping {
    std::vector<std::vector<std::string>> atom_site_ids;
    std::vector<std::string> model_ids;
};

[[nodiscard]] auto selected_structure_mapping(const ::gemmi::Structure& structure,
                                              const RecordSelection selection,
                                              const core::Molecule& molecule) -> BlockMapping {
    if (structure.models.size() != molecule.conformer_count()) {
        throw std::runtime_error{"structural model count does not match calculated conformers"};
    }

    BlockMapping mapping;
    mapping.atom_site_ids.reserve(structure.models.size());
    mapping.model_ids.reserve(structure.models.size());
    for (const auto& model : structure.models) {
        const auto selected = selection::SelectedModel{model, selection};
        if (selected.atoms().size() != molecule.atom_count()) {
            throw std::runtime_error{"structural atom count does not match calculated atoms"};
        }
        mapping.model_ids.push_back(std::to_string(model.num));
        auto& ids = mapping.atom_site_ids.emplace_back();
        ids.reserve(selected.atoms().size());
        for (std::size_t atom_index = 0; atom_index < selected.atoms().size(); ++atom_index) {
            const auto& source = *selected.atoms()[atom_index];
            const auto& atom = molecule.atom(atom_index);
            if (source.element.atomic_number() != atom.atomic_number() ||
                source.charge != atom.formal_charge() || source.name != atom.name()) {
                throw std::runtime_error{
                    "structural atom sequence does not match calculated atoms"};
            }
            ids.push_back(std::to_string(source.serial));
        }
    }
    return mapping;
}

[[nodiscard]] auto quote(const std::string_view value) -> std::string {
    return ::gemmi::cif::quote(std::string{value});
}

[[nodiscard]] auto block_name(const ImportedMoleculeRecord& record, const std::size_t index)
    -> std::string {
    auto candidate = record.identity.record_id.empty() ? std::string{record.molecule.name()}
                                                       : record.identity.record_id;
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
        return "SING";
    case core::BondOrder::DOUBLE:
        return "DOUB";
    case core::BondOrder::TRIPLE:
        return "TRIP";
    case core::BondOrder::AROMATIC:
        return "AROM";
    case core::BondOrder::UNKNOWN:
        throw std::runtime_error{"cannot write unknown bond order to generated mmCIF"};
    }
    throw std::runtime_error{"cannot write unsupported bond order to generated mmCIF"};
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

[[nodiscard]] auto atom_site_mapping(::gemmi::cif::Block& block, const core::Molecule& molecule,
                                     const RecordSelection record_selection) -> BlockMapping {
    const auto structure = ::gemmi::make_structure_from_block(block);
    auto mapping = selected_structure_mapping(structure, record_selection, molecule);
    auto table = block.find("_atom_site.", {"id", "type_symbol", "pdbx_PDB_model_num"});
    if (!table.ok()) {
        throw std::runtime_error{"mmCIF output block has no usable _atom_site category"};
    }

    BlockMapping row_mapping;
    std::string current_model;
    for (auto row : table) {
        const auto model = ::gemmi::cif::is_null(row[2]) ? "1" : ::gemmi::cif::as_string(row[2]);
        if (row_mapping.model_ids.empty() || current_model != model) {
            current_model = model;
            row_mapping.model_ids.push_back(model);
            row_mapping.atom_site_ids.emplace_back();
        }
        row_mapping.atom_site_ids.back().push_back(::gemmi::cif::as_string(row[0]));
    }

    if (row_mapping.model_ids.size() != molecule.conformer_count()) {
        throw std::runtime_error{"mmCIF model count does not match calculated conformers"};
    }
    for (const auto& ids : row_mapping.atom_site_ids) {
        if (ids.size() != molecule.atom_count()) {
            throw std::runtime_error{"mmCIF atom-site count does not match calculated atoms"};
        }
    }
    mapping.atom_site_ids = std::move(row_mapping.atom_site_ids);
    mapping.model_ids = std::move(row_mapping.model_ids);
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

[[nodiscard]] auto next_assignment_id(::gemmi::cif::Block& block) -> std::size_t {
    if (!block.has_mmcif_category(metadata_category)) {
        return 1;
    }
    auto table = block.find(metadata_category, {"id"});
    std::size_t next = 1;
    for (auto row : table) {
        try {
            next = std::max(
                next, static_cast<std::size_t>(std::stoull(::gemmi::cif::as_string(row[0]))) + 1);
        } catch (const std::exception&) {
            throw std::runtime_error{"existing mmCIF charge assignment ID is not numeric"};
        }
    }
    return next;
}

auto write_charges(::gemmi::cif::Block& block, const BlockMapping& mapping,
                   const core::Molecule& molecule,
                   const std::span<const charges::ChargeAssignment> assignments,
                   const charges::ChargeSet& charge_set, const std::string_view generator_name,
                   const std::string_view generator_version, const WriteMode mode) -> void {
    if (mode == WriteMode::replace) {
        erase_category(block, metadata_category);
        erase_category(block, charges_category);
    }
    ensure_dictionary(block);
    auto assignment_id = next_assignment_id(block);
    auto metadata =
        block.find_or_add(metadata_category, {"id", "type", "method", "parameter_set_id", "scope",
                                              "model_id", "software", "software_version"});
    auto charge_rows = block.find_or_add(charges_category, {"type_id", "atom_id", "charge"});
    metadata.ensure_loop();
    charge_rows.ensure_loop();

    for (const auto& assignment : assignments) {
        if (assignment.charges.size() != molecule.atom_count()) {
            throw std::runtime_error{"charge assignment size does not match mmCIF molecule"};
        }
        const auto conformer_index = assignment.target.conformer_index;
        const auto mapping_index = conformer_index.value_or(0);
        if (mapping_index >= mapping.atom_site_ids.size()) {
            throw std::runtime_error{"charge assignment conformer is missing from mmCIF block"};
        }
        const auto id = std::to_string(assignment_id++);
        const auto parameter = charge_set.parameter_set_id();
        metadata.append_row(
            {id, "empirical", quote(charge_set.method_id()),
             parameter.has_value() ? quote(*parameter) : ".",
             conformer_index.has_value() ? "model" : "topology",
             conformer_index.has_value() ? quote(mapping.model_ids[mapping_index]) : ".",
             quote(generator_name), generator_version.empty() ? "." : quote(generator_version)});
        for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
            charge_rows.append_row({id, quote(mapping.atom_site_ids[mapping_index][atom_index]),
                                    std::format("{:.4f}", assignment.charges[atom_index])});
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
                                 std::to_string(molecule_index)};
    }
    return result;
}

[[nodiscard]] auto selected_pdb_structure(const PdbSource& source) -> ::gemmi::Structure {
    auto structure = source.structure;
    for (auto& model : structure.models) {
        const auto selected = selection::SelectedModel{model, source.selection};
        std::unordered_set<const ::gemmi::Atom*> retained{selected.atoms().begin(),
                                                          selected.atoms().end()};
        for (auto& chain : model.chains) {
            for (auto& residue : chain.residues) {
                std::erase_if(residue.atoms, [&retained](const ::gemmi::Atom& atom) -> bool {
                    return !retained.contains(std::addressof(atom));
                });
            }
        }
        ::gemmi::remove_empty_children(model);
    }
    ::gemmi::setup_entities(structure);
    return structure;
}

} // namespace

MmcifWriter::MmcifWriter(std::ostream& output) : output_{std::addressof(output)} {}

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
                      generator_version, WriteMode::replace);
    }

    ::gemmi::cif::write_cif_to_stream(*output_, document);
    if (!*output_) {
        throw std::runtime_error{"failed to write mmCIF output"};
    }
}

auto MmcifWriter::write_pdb(const ImportedMoleculeRecord& record,
                            const charges::ChargeSet& charge_set, const PdbSource& source,
                            const std::string_view generator_name,
                            const std::string_view generator_version) const -> void {
    auto document = ::gemmi::make_mmcif_document(selected_pdb_structure(source));
    if (document.blocks.size() != 1) {
        throw std::runtime_error{"PDB conversion did not produce one mmCIF block"};
    }
    const auto assignments = assignments_for(charge_set, 0);
    auto& block = document.blocks.front();
    const auto mapping = atom_site_mapping(block, record.molecule, source.selection);
    write_charges(block, mapping, record.molecule, assignments, charge_set, generator_name,
                  generator_version, WriteMode::replace);

    ::gemmi::cif::write_cif_to_stream(*output_, document);
    if (!*output_) {
        throw std::runtime_error{"failed to write mmCIF output"};
    }
}

auto MmcifWriter::write_mmcif(const std::span<const ImportedMoleculeRecord> records,
                              const charges::ChargeSet& charge_set, const MmcifSource& source,
                              const std::string_view generator_name,
                              const std::string_view generator_version, const WriteMode mode) const
    -> void {
    if (source.document == nullptr || records.size() != source.block_indices.size()) {
        throw std::invalid_argument{"mmCIF source does not match imported records"};
    }
    auto document = *source.document;
    for (std::size_t record_index = 0; record_index < records.size(); ++record_index) {
        const auto block_index = source.block_indices[record_index];
        if (block_index >= document.blocks.size()) {
            throw std::runtime_error{"mmCIF source block index is out of range"};
        }
        const auto& record = records[record_index];
        const auto assignments = assignments_for(charge_set, record_index);
        auto& block = document.blocks[block_index];
        const auto mapping = atom_site_mapping(block, record.molecule, source.selection);
        write_charges(block, mapping, record.molecule, assignments, charge_set, generator_name,
                      generator_version, mode);
    }

    ::gemmi::cif::write_cif_to_stream(*output_, document);
    if (!*output_) {
        throw std::runtime_error{"failed to write mmCIF output"};
    }
}

} // namespace chargefw::adapters::gemmi::mmcif_output
