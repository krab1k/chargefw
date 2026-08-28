#include "bindings.h"

#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/methods/method_registry.h>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

namespace nb = nanobind;

namespace chargefw::python {
namespace {

auto option_value(const methods::MethodOptionValue& value) -> nb::object {
    return std::visit([](const auto& item) { return nb::cast(item); }, value);
}

auto method_descriptors() -> nb::list {
    auto result = nb::list{};
    for (const auto& method : methods::method_registry().methods()) {
        const auto& metadata = method->metadata();
        const auto requirements = method->requirements();
        auto descriptor = nb::dict{};
        descriptor["id"] = std::string{metadata.id};
        descriptor["name"] = std::string{metadata.name};
        descriptor["full_name"] = std::string{metadata.full_name};
        descriptor["publication"] = metadata.publication.has_value()
                                        ? nb::cast(std::string{*metadata.publication})
                                        : nb::none();
        descriptor["priority"] = metadata.priority;
        descriptor["requires_coordinates"] = requirements.coordinates;
        descriptor["supports_cutoff"] = requirements.resources.supports_cutoff;
        descriptor["supports_cover"] = requirements.resources.supports_cover;
        auto options = nb::list{};
        for (const auto& spec : method->option_schema()) {
            auto option = nb::dict{};
            option["id"] = std::string{spec.id};
            option["description"] = std::string{spec.description};
            option["type"] = nb::cast(spec.type);
            option["default"] = option_value(spec.default_value);
            auto choices = nb::list{};
            for (const auto& choice : spec.choices) {
                choices.append(option_value(choice));
            }
            option["choices"] = std::move(choices);
            option["minimum"] = spec.minimum.has_value() ? option_value(*spec.minimum) : nb::none();
            option["minimum_inclusive"] = spec.minimum_inclusive;
            option["maximum"] = spec.maximum.has_value() ? option_value(*spec.maximum) : nb::none();
            option["maximum_inclusive"] = spec.maximum_inclusive;
            options.append(std::move(option));
        }
        descriptor["options"] = std::move(options);
        result.append(std::move(descriptor));
    }
    return result;
}

} // namespace

void bind_methods(nb::module_& module) {
    nb::enum_<methods::MethodOptionType>(module, "MethodOptionType")
        .value("BOOLEAN", methods::MethodOptionType::boolean)
        .value("INTEGER", methods::MethodOptionType::integer)
        .value("FLOATING_POINT", methods::MethodOptionType::floating_point)
        .value("STRING", methods::MethodOptionType::string);
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
    module.def("_method_descriptors", &method_descriptors);
}

} // namespace chargefw::python
