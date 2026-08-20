#include <chargefw/methods/method.h>

#include <exception>
#include <string>

namespace chargefw::methods {
namespace {

auto add_invalid_options_issue(PrerequisiteResult& result, const std::exception& error) -> void {
    result.add(
        PrerequisiteIssue{.kind = PrerequisiteIssueKind::invalid_options, .message = error.what()});
}

} // namespace

auto Method::requires_parameters() const -> bool {
    return requirements().requires_parameters();
}

auto Method::check_method_prerequisites(const MethodPrerequisiteInput& input) const
    -> PrerequisiteResult {
    PrerequisiteResult result;

    try {
        validate_method_options(option_schema(), input.method_options);
    } catch (const std::exception& error) {
        add_invalid_options_issue(result, error);
    }

    const auto method_requirements = requirements();
    const auto& molecule = input.prepared_molecule.molecule();

    if (method_requirements.coordinates && !molecule.has_coordinates()) {
        result.add(PrerequisiteIssue{.kind = PrerequisiteIssueKind::missing_feature,
                                     .message = "method '" + std::string{id()} +
                                                "' requires coordinates"});
    }

    add_method_specific_prerequisite_issues(input, result);

    return result;
}

auto Method::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& /* unused */,
                                                     PrerequisiteResult& /* unused */) const
    -> void {}

} // namespace chargefw::methods
