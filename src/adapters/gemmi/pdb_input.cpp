#include <chargefw/adapters/gemmi/pdb_input.h>

#include "bonds.h"
#include "selection.h"
#include "structure_import.h"

#include <gemmi/pdb.hpp>

#include <algorithm>
#include <istream>
#include <iterator>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::adapters::gemmi::pdb_input {
namespace {

[[nodiscard]] auto optional_text(const std::string_view value) -> std::optional<std::string> {
    return value.empty() ? std::nullopt : std::optional<std::string>{value};
}

[[nodiscard]] auto optional_character(const char value) -> std::optional<std::string> {
    return value == '\0' || value == ' ' ? std::nullopt
                                         : std::optional<std::string>{std::string{value}};
}

[[nodiscard]] auto field(const std::string_view line, const std::size_t offset,
                         const std::size_t width) -> std::optional<std::string> {
    if (offset >= line.size()) {
        return std::nullopt;
    }
    const auto value = line.substr(offset, std::min(width, line.size() - offset));
    const auto first = value.find_first_not_of(' ');
    if (first == std::string_view::npos) {
        return std::nullopt;
    }
    const auto last = value.find_last_not_of(' ');
    return std::string{value.substr(first, last - first + 1)};
}

struct PdbSourceModel {
    std::optional<std::string> id;
    std::vector<SourceAtomReference> sites;
};

struct PdbSource {
    std::vector<PdbSourceModel> models;
    SourceConnectivity connectivity = SourceConnectivity::absent;
};

[[nodiscard]] auto parse_source(const std::string& contents) -> PdbSource {
    auto result = PdbSource{};
    auto* current = static_cast<PdbSourceModel*>(nullptr);
    auto source_position = std::size_t{0};
    auto input = std::istringstream{contents};
    auto line = std::string{};
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        const auto view = std::string_view{line};
        if (view.starts_with("CONECT") || view.starts_with("LINK  ") ||
            view.starts_with("SSBOND")) {
            result.connectivity = SourceConnectivity::present;
        }
        if (view.starts_with("MODEL ")) {
            result.models.push_back(PdbSourceModel{.id = field(view, 10, 4), .sites = {}});
            current = std::addressof(result.models.back());
            continue;
        }
        if (view.starts_with("ENDMDL")) {
            current = nullptr;
            continue;
        }
        if (!view.starts_with("ATOM  ") && !view.starts_with("HETATM")) {
            continue;
        }
        if (current == nullptr) {
            result.models.push_back(PdbSourceModel{.id = std::nullopt, .sites = {}});
            current = std::addressof(result.models.back());
        }
        current->sites.push_back(
            SourceAtomReference{.position = source_position++,
                                .id = field(view, 6, 5),
                                .structural_labels = SourceStructuralLabels{
                                    .author = SourceHierarchyLabels{.atom = field(view, 12, 4),
                                                                    .residue = field(view, 17, 3),
                                                                    .chain = field(view, 21, 1),
                                                                    .sequence = field(view, 22, 4)},
                                    .label = {},
                                    .entity = std::nullopt,
                                    .insertion_code = field(view, 26, 1),
                                    .alternate_location = field(view, 16, 1),
                                    .segment = field(view, 72, 4)}});
    }
    return result;
}

[[nodiscard]] auto selected_labels(const selection::SelectedSite& site) -> SourceStructuralLabels {
    return SourceStructuralLabels{
        .author = SourceHierarchyLabels{.atom = optional_text(site.atom->name),
                                        .residue = optional_text(site.residue->name),
                                        .chain = optional_text(site.chain_name),
                                        .sequence = optional_text(site.residue->seqid.num.str())},
        .label = {},
        .entity = std::nullopt,
        .insertion_code = optional_character(site.residue->seqid.icode),
        .alternate_location = optional_character(site.atom->altloc),
        .segment = optional_text(site.residue->segment)};
}

[[nodiscard]] auto source_models(const ::gemmi::Structure& structure, const InputOptions options,
                                 const std::vector<PdbSourceModel>& source)
    -> std::vector<structure_import::SourceModelMapping> {
    const auto retained_count =
        options.conformers == ConformerSelection::all ? structure.models.size() : std::size_t{1};
    if (source.size() != structure.models.size()) {
        throw std::runtime_error{"PDB source mapping does not match parsed models"};
    }
    auto result = std::vector<structure_import::SourceModelMapping>{};
    result.reserve(retained_count);
    for (std::size_t model_index = 0; model_index < retained_count; ++model_index) {
        const auto& model = structure.models[model_index];
        const auto selected = selection::SelectedModel{model, options.selection};
        auto sites = std::vector<SourceAtomReference>{};
        sites.reserve(selected.sites().size());
        auto used = std::vector<bool>(source[model_index].sites.size(), false);
        for (const auto& site : selected.sites()) {
            const auto labels = selected_labels(site);
            auto match = std::optional<std::size_t>{};
            for (std::size_t index = 0; index < source[model_index].sites.size(); ++index) {
                if (!used[index] && source[model_index].sites[index].structural_labels == labels) {
                    if (match.has_value()) {
                        throw std::runtime_error{
                            "PDB source contains ambiguous selected atom identity"};
                    }
                    match = index;
                }
            }
            if (!match.has_value()) {
                throw std::runtime_error{"PDB source mapping does not match selected atoms"};
            }
            used[*match] = true;
            sites.push_back(source[model_index].sites[*match]);
        }
        result.push_back(structure_import::SourceModelMapping{
            .conformer = SourceConformerReference{
                .position = model_index, .id = source[model_index].id, .sites = std::move(sites)}});
    }
    return result;
}

} // namespace

PdbReader::PdbReader(std::istream& input, std::string source,
                     const ::chargefw::adapters::gemmi::InputOptions options)
    : options_{options} {
    std::string contents{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    if (input.bad()) {
        throw std::runtime_error{"failed to read PDB input"};
    }

    structure_ = ::gemmi::read_pdb_string(contents, source);
    if (structure_.models.empty()) {
        throw std::runtime_error{"structural input contains no models"};
    }
    const auto name = structure_.name.empty() ? source : structure_.name;
    const auto selected = selection::SelectedModel{structure_.models.front(), options.selection};
    auto explicit_bonds = bonds::explicit_pdb(structure_, selected);
    const auto parsed_source = parse_source(contents);
    auto mappings = source_models(structure_, options, parsed_source.models);
    record_ = structure_import::make_record(
        structure_,
        MoleculeRecordIdentity{.source = std::move(source), .record_index = 0, .record_id = name},
        options.selection, options.bond_strategy, options.conformers, std::move(explicit_bonds),
        name, MolecularSourceFormat::pdb, std::move(mappings), parsed_source.connectivity);
}

auto PdbReader::next() -> std::optional<ImportedMoleculeRecord> {
    return std::exchange(record_, std::nullopt);
}

auto PdbReader::source_structure() const -> const ::gemmi::Structure& {
    return structure_;
}

auto PdbReader::options() const noexcept -> ::chargefw::adapters::gemmi::InputOptions {
    return options_;
}

} // namespace chargefw::adapters::gemmi::pdb_input
