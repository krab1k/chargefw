#pragma once

#include <chargefw/adapters/molecule_record.h>

#include <cstddef>
#include <optional>
#include <string>

#include <gemmi/cifdoc.hpp>

namespace chargefw::adapters::gemmi::mmcif_labels {

struct SourceLabelColumns {
    std::size_t label_atom;
    std::size_t label_residue;
    std::size_t label_chain;
    std::size_t label_sequence;
    std::size_t author_atom;
    std::size_t author_residue;
    std::size_t author_chain;
    std::size_t author_sequence;
    std::size_t insertion_code;
    std::size_t alternate_location;
    std::size_t entity;
};

[[nodiscard]] inline auto source_value(const ::gemmi::cif::Table::Row& row, const std::size_t index)
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

[[nodiscard]] inline auto decode_source_labels(const ::gemmi::cif::Table::Row& row,
                                               const SourceLabelColumns columns)
    -> SourceStructuralLabels {
    return SourceStructuralLabels{
        .author = SourceHierarchyLabels{.atom = source_value(row, columns.author_atom),
                                        .residue = source_value(row, columns.author_residue),
                                        .chain = source_value(row, columns.author_chain),
                                        .sequence = source_value(row, columns.author_sequence)},
        .label = SourceHierarchyLabels{.atom = source_value(row, columns.label_atom),
                                       .residue = source_value(row, columns.label_residue),
                                       .chain = source_value(row, columns.label_chain),
                                       .sequence = source_value(row, columns.label_sequence)},
        .entity = source_value(row, columns.entity),
        .insertion_code = source_value(row, columns.insertion_code),
        .alternate_location = source_value(row, columns.alternate_location),
        .segment = std::nullopt};
}

} // namespace chargefw::adapters::gemmi::mmcif_labels
