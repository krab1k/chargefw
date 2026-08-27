#include "bindings.h"

#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_prerequisites.h>

#include <nanobind/nanobind.h>

namespace nb = nanobind;

namespace chargefw::python {

void bind_methods(nb::module_& module) {
    nb::enum_<methods::PrerequisiteIssueKind>(module, "PrerequisiteIssueKind")
        .value("INVALID_OPTIONS", methods::PrerequisiteIssueKind::invalid_options)
        .value("MISSING_FEATURE", methods::PrerequisiteIssueKind::missing_feature)
        .value("INVALID_GEOMETRY", methods::PrerequisiteIssueKind::invalid_geometry)
        .value("UNSUPPORTED_MOLECULE", methods::PrerequisiteIssueKind::unsupported_molecule)
        .value("MISSING_PARAMETERS", methods::PrerequisiteIssueKind::missing_parameters)
        .value("PARAMETER_CLASSIFICATION_FAILED",
               methods::PrerequisiteIssueKind::parameter_classification_failed);
    nb::enum_<methods::ExecutionAvailability>(module, "ExecutionAvailability")
        .value("AVAILABLE", methods::ExecutionAvailability::available)
        .value("AVAILABLE_WITH_WARNING", methods::ExecutionAvailability::available_with_warning)
        .value("UNSUPPORTED", methods::ExecutionAvailability::unsupported);
    nb::enum_<methods::ExecutionIssueKind>(module, "ExecutionIssueKind")
        .value("RESOURCE_THRESHOLD_EXCEEDED",
               methods::ExecutionIssueKind::resource_threshold_exceeded)
        .value("UNSUPPORTED_EXECUTION_MODE",
               methods::ExecutionIssueKind::unsupported_execution_mode);
}

} // namespace chargefw::python
