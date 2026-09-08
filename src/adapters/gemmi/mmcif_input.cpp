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
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace chargefw::adapters::gemmi::mmcif_input {
namespace {

auto validate_atom_site_ids(::gemmi::cif::Block& block) -> void {
    auto atom_sites = block.find("_atom_site.", {"id"});
    std::unordered_set<int> ids;
    ids.reserve(atom_sites.length());
    for (auto row : atom_sites) {
        const auto source_id = ::gemmi::cif::as_string(row[0]);
        int id = 0;
        const auto [end, error] =
            std::from_chars(source_id.data(), source_id.data() + source_id.size(), id);
        if (error != std::errc{} || end != source_id.data() + source_id.size() ||
            std::to_string(id) != source_id) {
            throw std::runtime_error{"mmCIF _atom_site.id must be a unique canonical integer; "
                                     "unsupported value '" +
                                     source_id + "'"};
        }
        if (!ids.insert(id).second) {
            throw std::runtime_error{"mmCIF _atom_site.id must be a unique canonical integer; "
                                     "duplicate value '" +
                                     source_id + "'"};
        }
    }
}

} // namespace

MmcifReader::MmcifReader(std::istream& input, std::string source,
                         const ::chargefw::adapters::gemmi::InputOptions options)
    : source_{std::move(source)}, options_{options} {
    std::string contents{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    if (input.bad()) {
        throw std::runtime_error{"failed to read mmCIF input"};
    }

    document_ = std::make_shared<::gemmi::cif::Document>(
        ::gemmi::cif::read_memory(contents.data(), contents.size(), source_.c_str()));

    const auto has_coordinates = [](const ::gemmi::cif::Block& block) -> bool {
        return block.has_tag("_atom_site.id");
    };
    if (std::ranges::none_of(document_->blocks, has_coordinates)) {
        throw std::runtime_error{"mmCIF input contains no coordinate data blocks"};
    }
}

auto MmcifReader::next() -> std::optional<ImportedMoleculeRecord> {
    while (block_index_ < document_->blocks.size()) {
        const auto current_block_index = block_index_;
        auto& block = document_->blocks[block_index_++];
        if (!block.has_tag("_atom_site.id")) {
            continue;
        }

        validate_atom_site_ids(block);
        const auto structure = ::gemmi::make_structure_from_block(block);
        if (structure.models.empty()) {
            throw std::runtime_error{"structural input contains no models"};
        }
        const auto selected =
            selection::SelectedModel{structure.models.front(), options_.selection};
        auto explicit_bonds = bonds::explicit_mmcif(structure, block, selected);
        const auto current_record_index = record_index_++;
        auto record = structure_import::make_record(
            structure,
            MoleculeRecordIdentity{
                .source = source_, .record_index = current_record_index, .record_id = block.name},
            options_.selection, options_.bond_strategy, options_.conformers,
            std::move(explicit_bonds), structure.name.empty() ? block.name : structure.name);
        source_block_indices_.push_back(current_block_index);
        return record;
    }

    return std::nullopt;
}

auto MmcifReader::source_document() const -> std::shared_ptr<const ::gemmi::cif::Document> {
    return document_;
}

auto MmcifReader::source_block_indices() const -> const std::vector<std::size_t>& {
    return source_block_indices_;
}

auto MmcifReader::options() const noexcept -> ::chargefw::adapters::gemmi::InputOptions {
    return options_;
}

} // namespace chargefw::adapters::gemmi::mmcif_input
