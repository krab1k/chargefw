#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>

#include <chargefw/adapters/molecule_record.h>

namespace chargefw::adapters::detail {

[[nodiscard]] inline auto usable_source_name(const std::optional<std::string>& name) -> bool {
    return name && !name->empty() && *name != "." && *name != "?";
}

struct ComponentAtom {
    std::string name;
    std::size_t index;
};

[[nodiscard]] inline auto
canonical_component_name(const std::optional<SourceStructuralLabels>& labels,
                         const std::string_view fallback) -> std::string_view {
    if (labels && usable_source_name(labels->label.residue)) {
        return *labels->label.residue;
    }
    if (labels && usable_source_name(labels->author.residue)) {
        return *labels->author.residue;
    }
    return fallback;
}

[[nodiscard]] inline auto canonical_atom_name(const std::optional<SourceStructuralLabels>& labels,
                                              const std::string_view fallback) -> std::string_view {
    if (labels && usable_source_name(labels->label.atom)) {
        return *labels->label.atom;
    }
    if (labels && usable_source_name(labels->author.atom)) {
        return *labels->author.atom;
    }
    return fallback;
}

struct ComponentView {
    std::string_view name;
    std::span<const ComponentAtom> atoms;

    [[nodiscard]] auto find_atom(const std::string_view atom_name) const
        -> std::optional<std::size_t> {
        for (const auto& atom : atoms) {
            if (atom.name == atom_name) {
                return atom.index;
            }
        }
        return std::nullopt;
    }
};

} // namespace chargefw::adapters::detail
