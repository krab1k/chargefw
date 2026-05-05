#include <chargefw/core/periodic_table.h>
#include <chargefw/methods/method.h>

#include <exception>
#include <string>

namespace chargefw::methods {
namespace {

auto add_invalid_options_issue(PrerequisiteResult& result, const std::exception& error) -> void {
    result.add(
        PrerequisiteIssue{.kind = PrerequisiteIssueKind::invalid_options, .message = error.what()});
}

auto add_missing_element_properties_issues(const Method& method, const core::Molecule& molecule,
                                           PrerequisiteResult& result) -> void {
    const auto& table = core::periodic_table();

    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto& atom = molecule.atom(atom_index);

        if (!table.contains(atom.atomic_number())) {
            result.add(PrerequisiteIssue{
                .kind = PrerequisiteIssueKind::unsupported_molecule,
                .message = "method '" + std::string{method.id()} +
                           "' requires element properties for atom " + std::to_string(atom_index) +
                           " with atomic number " + std::to_string(atom.atomic_number()),
                .atom_index = atom_index});
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

    if (method_requirements.element_properties) {
        add_missing_element_properties_issues(*this, molecule, result);
    }

    if (method_requirements.resources.reject_large_without_reduction &&
        molecule.atom_count() > method_requirements.resources.large_molecule_atom_threshold) {
        result.add(PrerequisiteIssue{.kind = PrerequisiteIssueKind::resource_limit,
                                     .message = "method '" + std::string{id()} +
                                                "' is not suitable for molecule '" +
                                                std::string{molecule.name()} + "' with " +
                                                std::to_string(molecule.atom_count()) + " atoms"});
    }

    add_method_specific_prerequisite_issues(input, result);

    return result;
}

auto Method::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& /* unused */,
                                                     PrerequisiteResult& /* unused */) const
    -> void {}

} // namespace chargefw::methods