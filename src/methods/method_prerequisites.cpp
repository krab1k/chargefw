#include <chargefw/methods/method.h>
#include <chargefw/methods/method_prerequisites.h>

#include "core/diagnostic_description.h"

#include <string>
#include <utility>

namespace chargefw::methods {
namespace {

[[nodiscard]] auto with_molecule_prefix(const core::Molecule& molecule,
                                        const std::size_t molecule_index,
                                        const PrerequisiteIssue& issue) -> PrerequisiteIssue {
    return PrerequisiteIssue{.kind = issue.kind,
                             .message =
                                 core::detail::molecule_description(molecule, molecule_index) +
                                 ": " + issue.message,
                             .molecule_index = molecule_index,
                             .atom_index = issue.atom_index,
                             .bond_index = issue.bond_index,
                             .conformer_index = issue.conformer_index};
}

} // namespace

auto PrerequisiteResult::ok() const noexcept -> bool {
    return issues_.empty();
}

PrerequisiteResult::operator bool() const noexcept {
    return ok();
}

auto PrerequisiteResult::issues() const noexcept -> std::span<const PrerequisiteIssue> {
    return issues_;
}

auto PrerequisiteResult::add(PrerequisiteIssue issue) -> void {
    issues_.push_back(std::move(issue));
}

auto check_method_prerequisites(const Method& method,
                                const features::PreparedMoleculeCollection& molecules,
                                const MethodOptions& method_options) -> PrerequisiteResult {
    auto result = PrerequisiteResult{};

    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        const auto molecule_result = method.check_method_prerequisites(
            {.prepared_molecule = molecules[molecule_index], .method_options = method_options});

        for (const auto& issue : molecule_result.issues()) {
            result.add(
                with_molecule_prefix(molecules[molecule_index].molecule(), molecule_index, issue));
        }
    }

    return result;
}

} // namespace chargefw::methods
