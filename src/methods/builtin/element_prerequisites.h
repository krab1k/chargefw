#pragma once

#include "core/diagnostic_description.h"

#include <chargefw/core/periodic_table.h>
#include <chargefw/methods/method_prerequisites.h>

#include <concepts>
#include <cstddef>
#include <string>
#include <string_view>

namespace chargefw::methods::builtin::detail {

template <typename ElementSupported>
    requires std::predicate<ElementSupported, const core::Element&>

auto add_element_prerequisite_issues(const MethodPrerequisiteInput& input,
                                     PrerequisiteResult& result, const std::string_view requirement,
                                     const ElementSupported& element_supported) -> void {
    const auto& table = core::periodic_table();
    const auto& molecule = input.prepared_molecule.molecule();

    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto atomic_number = molecule.atom(atom_index).atomic_number();

        if (!table.contains(atomic_number) || !element_supported(table.element(atomic_number))) {
            result.add(
                PrerequisiteIssue{.kind = PrerequisiteIssueKind::unsupported_molecule,
                                  .message = std::string{requirement} + ": " +
                                             core::detail::atom_description(molecule, atom_index),
                                  .atom_index = atom_index});
        }
    }
}

} // namespace chargefw::methods::builtin::detail
