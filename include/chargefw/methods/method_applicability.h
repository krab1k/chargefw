#pragma once

#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/parameters/classification/classification_result.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <optional>
#include <span>
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
};

struct ApplicableMethod {
    const Method* method = nullptr;
    const parameters::ParameterSet* parameter_set = nullptr;

    MethodOptions method_options;
    std::vector<parameters::ParameterClassification> classifications;

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
