#pragma once

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/parameters/classification/classification_result.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace chargefw::methods {

class Method;

struct ApplicabilityRequest {
    // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
    const features::PreparedMoleculeCollection& molecules;
    std::span<const Method* const> methods;
    std::span<const parameters::ParameterSet> parameter_sets;
    // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
    parameters::ClassificationOptions classification_options{};
    calculation::ResourcePolicy resource_policy{};
};

enum class ExecutionAvailability : std::uint8_t {
    available,
    available_with_warning,
    unsupported,
};

enum class ExecutionIssueKind : std::uint8_t {
    resource_threshold_exceeded,
    unsupported_execution_mode,
};

struct ExecutionIssue {
    ExecutionIssueKind kind;
    std::string message;
    std::optional<std::size_t> molecule_index = std::nullopt;
};

struct ExecutionAssessment {
    calculation::ExecutionMode mode;
    ExecutionAvailability availability;
    std::vector<ExecutionIssue> issues;
};

struct ApplicableMethod {
    const Method* method = nullptr;
    const parameters::ParameterSet* parameter_set = nullptr;

    MethodOptions method_options;
    std::vector<parameters::ParameterClassification> classifications;
    std::vector<ExecutionAssessment> execution_assessments;

    [[nodiscard]] auto uses_parameters() const noexcept -> bool {
        return parameter_set != nullptr;
    }
};

struct RejectedCandidate {
    std::size_t method_index = 0;
    std::optional<std::size_t> parameter_set_index;
    std::vector<PrerequisiteIssue> issues;
};

struct ApplicabilityResult {
    std::vector<ApplicableMethod> applicable;
    std::vector<RejectedCandidate> rejected;

    [[nodiscard]] auto empty() const noexcept -> bool {
        return applicable.empty();
    }
};

[[nodiscard]] auto find_applicable_methods(const ApplicabilityRequest& request)
    -> ApplicabilityResult;

} // namespace chargefw::methods
