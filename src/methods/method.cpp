#include <chargefw/methods/method.h>

#include <chargefw/features/conformer_features.h>

#include <exception>
#include <string>

namespace chargefw::methods {
namespace {

auto add_invalid_options_issue(PrerequisiteResult& result, const std::exception& error) -> void {
    result.add(
        PrerequisiteIssue{.kind = PrerequisiteIssueKind::invalid_options, .message = error.what()});
}

auto add_coordinate_prerequisite_issues(const Method& method, const core::Molecule& molecule,
                                        PrerequisiteResult& result) -> void {
    for (std::size_t conformer_index = 0; conformer_index < molecule.conformer_count();
         ++conformer_index) {
        const features::ConformerFeatures geometry{molecule, conformer_index};

        if (const auto nonfinite_atom = geometry.first_nonfinite_atom_index()) {
            result.add(PrerequisiteIssue{.kind = PrerequisiteIssueKind::invalid_geometry,
                                         .message = "method '" + std::string{method.id()} +
                                                    "' requires finite coordinates; conformer " +
                                                    std::to_string(conformer_index),
                                         .atom_index = *nonfinite_atom});
            continue;
        }

        if (const auto coincident = geometry.coincident_atom_indices()) {
            result.add(PrerequisiteIssue{
                .kind = PrerequisiteIssueKind::invalid_geometry,
                .message = "method '" + std::string{method.id()} +
                           "' requires distinct atom coordinates; conformer " +
                           std::to_string(conformer_index) + ", atoms " +
                           std::to_string(coincident->first) + " and " +
                           std::to_string(coincident->second) + " share coordinates",
                .atom_index = coincident->second});
        }
    }
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
    if (method_requirements.coordinates && molecule.has_coordinates()) {
        add_coordinate_prerequisite_issues(*this, molecule, result);
    }

    add_method_specific_prerequisite_issues(input, result);

    return result;
}

auto Method::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& /* unused */,
                                                     PrerequisiteResult& /* unused */) const
    -> void {}

} // namespace chargefw::methods
