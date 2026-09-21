#include <chargefw/adapters/gemmi/mmcif_input.h>

#include "bonds.h"
#include "selection.h"
#include "structure_import.h"

#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>

#include <algorithm>
#include <charconv>
#include <istream>
#include <iterator>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace chargefw::adapters::gemmi::mmcif_input {
namespace {

[[nodiscard]] auto parser_block(const ::gemmi::cif::Block& source) -> ::gemmi::cif::Block {
    auto result = source;
    auto atom_sites = result.find("_atom_site.", {"id"});
    if (atom_sites.length() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error{"mmCIF contains too many atom sites for Gemmi"};
    }
    auto ids = std::unordered_set<std::string>{};
    ids.reserve(atom_sites.length());
    auto position = std::size_t{0};
    for (auto row : atom_sites) {
        const auto source_id = ::gemmi::cif::as_string(row[0]);
        if (source_id.empty() || source_id == "." || source_id == "?") {
            throw std::runtime_error{"mmCIF _atom_site.id must not be missing or unknown"};
        }
        if (!ids.insert(source_id).second) {
            throw std::runtime_error{"mmCIF _atom_site.id must be unique; duplicate value '" +
                                     source_id + "'"};
        }
        row[0] = std::to_string(++position);
    }
    return result;
}

[[nodiscard]] auto source_value(const ::gemmi::cif::Table::Row& row, const std::size_t index)
    -> std::optional<std::string> {
    if (!row.has(index)) {
        return std::nullopt;
    }
    const auto& token = row[index];
    if (token == "." || token == "?") {
        return token;
    }
    return ::gemmi::cif::as_string(token);
}

struct MmcifSourceSite {
    int model_number = 1;
    std::optional<std::string> model_id;
    SourceAtomReference reference;
};

[[nodiscard]] auto source_sites(::gemmi::cif::Block& block) -> std::vector<MmcifSourceSite> {
    auto atom_sites =
        block.find("_atom_site.", {"id", "?label_atom_id", "?label_comp_id", "?label_asym_id",
                                   "?label_seq_id", "?auth_atom_id", "?auth_comp_id",
                                   "?auth_asym_id", "?auth_seq_id", "?pdbx_PDB_ins_code",
                                   "?label_alt_id", "?label_entity_id", "?pdbx_PDB_model_num"});
    auto result = std::vector<MmcifSourceSite>{};
    result.reserve(atom_sites.length());
    auto position = std::size_t{0};
    for (const auto row : atom_sites) {
        const auto id = ::gemmi::cif::as_string(row[0]);
        const auto model_id = source_value(row, 12);
        auto model_number = 1;
        if (model_id.has_value()) {
            const auto [end, error] = std::from_chars(
                model_id->data(), model_id->data() + model_id->size(), model_number);
            if (error != std::errc{} || end != model_id->data() + model_id->size()) {
                throw std::runtime_error{"mmCIF model IDs must be integers"};
            }
        }
        result.push_back(MmcifSourceSite{
            .model_number = model_number,
            .model_id = model_id,
            .reference = SourceAtomReference{
                .position = position++,
                .id = id,
                .structural_labels = SourceStructuralLabels{
                    .author = SourceHierarchyLabels{.atom = source_value(row, 5),
                                                    .residue = source_value(row, 6),
                                                    .chain = source_value(row, 7),
                                                    .sequence = source_value(row, 8)},
                    .label = SourceHierarchyLabels{.atom = source_value(row, 1),
                                                   .residue = source_value(row, 2),
                                                   .chain = source_value(row, 3),
                                                   .sequence = source_value(row, 4)},
                    .entity = source_value(row, 11),
                    .insertion_code = source_value(row, 9),
                    .alternate_location = source_value(row, 10),
                    .segment = std::nullopt}}});
    }
    return result;
}

[[nodiscard]] auto make_source_mappings(::gemmi::cif::Block& block,
                                        const ::gemmi::Structure& structure,
                                        const InputOptions options)
    -> std::vector<structure_import::SourceModelMapping> {
    const auto sites = source_sites(block);
    const auto retained_count =
        options.conformers == ConformerSelection::all ? structure.models.size() : std::size_t{1};
    auto result = std::vector<structure_import::SourceModelMapping>{};
    result.reserve(retained_count);
    for (std::size_t model_index = 0; model_index < retained_count; ++model_index) {
        const auto& model = structure.models[model_index];
        const auto selected = selection::SelectedModel{model, options.selection};
        auto mapped = std::vector<SourceAtomReference>{};
        mapped.reserve(selected.atoms().size());
        auto model_id = std::optional<std::string>{};
        auto model_id_initialized = false;
        for (const auto* atom : selected.atoms()) {
            if (atom->serial <= 0 || static_cast<std::size_t>(atom->serial) > sites.size()) {
                throw std::runtime_error{"mmCIF source mapping does not match selected atoms"};
            }
            const auto& source = sites[static_cast<std::size_t>(atom->serial) - 1];
            if (source.model_number != model.num) {
                throw std::runtime_error{"mmCIF source mapping does not match selected atoms"};
            }
            if (model_id_initialized && model_id != source.model_id) {
                throw std::runtime_error{"mmCIF model contains inconsistent source model IDs"};
            }
            model_id = source.model_id;
            model_id_initialized = true;
            mapped.push_back(source.reference);
        }
        result.push_back(structure_import::SourceModelMapping{
            .conformer = SourceConformerReference{
                .position = model_index, .id = std::move(model_id), .sites = std::move(mapped)}});
    }
    return result;
}

[[nodiscard]] auto source_connectivity(const ::gemmi::cif::Block& block) -> SourceConnectivity {
    return block.has_mmcif_category("_chem_comp_bond.") || block.has_mmcif_category("_struct_conn.")
               ? SourceConnectivity::present
               : SourceConnectivity::absent;
}

} // namespace

MmcifReader::MmcifReader(std::istream& input, std::string source,
                         const ::chargefw::adapters::gemmi::InputOptions options)
    : source_{std::move(source)}, options_{options} {
    std::string contents{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    if (input.bad()) {
        throw std::runtime_error{"failed to read mmCIF input"};
    }

    document_ = ::gemmi::cif::read_memory(contents.data(), contents.size(), source_.c_str());

    const auto has_coordinates = [](const ::gemmi::cif::Block& block) -> bool {
        return block.has_tag("_atom_site.id");
    };
    if (std::ranges::none_of(document_.blocks, has_coordinates)) {
        throw std::runtime_error{"mmCIF input contains no coordinate data blocks"};
    }
}

auto MmcifReader::next() -> std::optional<ImportedMoleculeRecord> {
    while (block_index_ < document_.blocks.size()) {
        auto& block = document_.blocks[block_index_++];
        if (!block.has_tag("_atom_site.id")) {
            continue;
        }

        auto normalized = parser_block(block);
        const auto structure = ::gemmi::make_structure_from_block(normalized);
        if (structure.models.empty()) {
            throw std::runtime_error{"structural input contains no models"};
        }
        const auto selected =
            selection::SelectedModel{structure.models.front(), options_.selection};
        auto explicit_bonds = bonds::explicit_mmcif(structure, block, selected);
        auto source_models = make_source_mappings(block, structure, options_);
        const auto current_record_index = record_index_++;
        auto record = structure_import::make_record(
            structure,
            MoleculeRecordIdentity{
                .source = source_, .record_index = current_record_index, .record_id = block.name},
            options_.selection, options_.bond_strategy, options_.conformers,
            std::move(explicit_bonds), structure.name.empty() ? block.name : structure.name,
            MolecularSourceFormat::mmcif, std::move(source_models), source_connectivity(block));
        return record;
    }

    return std::nullopt;
}

auto MmcifReader::options() const noexcept -> ::chargefw::adapters::gemmi::InputOptions {
    return options_;
}

} // namespace chargefw::adapters::gemmi::mmcif_input
