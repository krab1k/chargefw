#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/methods/method_options.h>

#include <optional>
#include <span>
#include <string>
#include <vector>

namespace chargefw::methods {

enum class PrerequisiteIssueKind {
    invalid_options,
    missing_feature,
    unsupported_molecule,
    resource_limit,
    missing_parameters,
    parameter_classification_failed,
};

struct PrerequisiteIssue {
    PrerequisiteIssueKind kind;
    std::string message;
    std::optional<std::size_t> atom_index = std::nullopt;
    std::optional<std::size_t> bond_index = std::nullopt;
};

class PrerequisiteResult {
  public:
    [[nodiscard]] auto ok() const noexcept -> bool;
    [[nodiscard]] explicit operator bool() const noexcept;

    [[nodiscard]] auto issues() const noexcept -> std::span<const PrerequisiteIssue>;

    auto add(PrerequisiteIssue issue) -> void;

  private:
    std::vector<PrerequisiteIssue> issues_;
};

struct MethodPrerequisiteInput {
    const core::Molecule& molecule;
    const MethodOptions& method_options;
};

} // namespace chargefw::methods