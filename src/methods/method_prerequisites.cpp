#include <chargefw/methods/method.h>
#include <chargefw/methods/method_prerequisites.h>

#include <string>
#include <utility>

namespace chargefw::methods {
namespace {

[[nodiscard]] auto with_molecule_prefix(std::size_t molecule_index, const PrerequisiteIssue& issue)
    -> PrerequisiteIssue {
    return PrerequisiteIssue{.kind = issue.kind,
                             .message = "molecule " + std::to_string(molecule_index) + ": " +
                                        issue.message,
                             .atom_index = issue.atom_index,
                             .bond_index = issue.bond_index};
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
            result.add(with_molecule_prefix(molecule_index, issue));
        }
    }

    return result;
}

} // namespace chargefw::methods