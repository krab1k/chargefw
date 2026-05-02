#pragma once

#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/parameters/classification/classification_result.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <optional>
#include <vector>

namespace chargefw::methods {

class Method;

struct ParameterPrerequisiteInput {
    // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
    const features::PreparedMolecule& prepared_molecule;
    const parameters::ParameterSet& parameter_set;
    // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
    parameters::ClassificationOptions classification_options = {};
};

struct ParameterPrerequisiteResult {
    std::vector<PrerequisiteIssue> issues;
    std::optional<parameters::ParameterClassification> classification;

    [[nodiscard]] auto ok() const noexcept -> bool {
        return issues.empty();
    }

    [[nodiscard]] explicit operator bool() const noexcept {
        return ok();
    }
};

struct CollectionParameterPrerequisiteResult {
    std::vector<PrerequisiteIssue> issues;
    std::vector<parameters::ParameterClassification> classifications;

    [[nodiscard]] auto ok() const noexcept -> bool {
        return issues.empty();
    }

    [[nodiscard]] explicit operator bool() const noexcept {
        return ok();
    }
};

[[nodiscard]] auto check_parameter_prerequisites(const Method& method,
                                                 const ParameterPrerequisiteInput& input)
    -> ParameterPrerequisiteResult;

[[nodiscard]] auto
check_parameter_prerequisites(const Method& method,
                              const features::PreparedMoleculeCollection& molecules,
                              const parameters::ParameterSet& parameter_set,
                              const parameters::ClassificationOptions& classification_options = {})
    -> CollectionParameterPrerequisiteResult;

} // namespace chargefw::methods