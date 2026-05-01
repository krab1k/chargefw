#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/features/topology_features.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/parameters/classification_result.h>
#include <chargefw/parameters/parameter_classification.h>
#include <chargefw/parameters/parameter_set.h>

#include <optional>
#include <vector>

namespace chargefw::methods {

struct ParameterPrerequisiteInput {
    const core::Molecule& molecule;
    const features::TopologyFeatures& topology;
    const parameters::ParameterSet& parameter_set;
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

[[nodiscard]] auto check_parameter_prerequisites(const Method& method,
                                                 const ParameterPrerequisiteInput& input)
    -> ParameterPrerequisiteResult;

} // namespace chargefw::methods