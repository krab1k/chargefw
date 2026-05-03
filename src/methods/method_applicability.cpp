#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/parameter_prerequisites.h>

#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::methods {
namespace {

[[nodiscard]] auto make_issue(const PrerequisiteIssueKind kind, std::string message)
    -> PrerequisiteIssue {
    return PrerequisiteIssue{.kind = kind, .message = std::move(message)};
}

auto add_rejected(ApplicabilityResult& result, const std::size_t method_index,
                  const std::optional<std::size_t> parameter_set_index,
                  std::vector<PrerequisiteIssue> issues) -> void {
    result.rejected.push_back(RejectedCandidate{.method_index = method_index,
                                                .parameter_set_index = parameter_set_index,
                                                .issues = std::move(issues)});
}

[[nodiscard]] auto copy_issues(std::span<const PrerequisiteIssue> issues)
    -> std::vector<PrerequisiteIssue> {
    return {issues.begin(), issues.end()};
}

[[nodiscard]] auto has_parameter_classification(const MethodRequirements& requirements) noexcept
    -> bool {
    return requirements.requires_atom_parameters() || requirements.requires_bond_parameters();
}

[[nodiscard]] auto
valid_classification_count(const MethodRequirements& requirements,
                           const features::PreparedMoleculeCollection& molecules,
                           const CollectionParameterPrerequisiteResult& result) noexcept -> bool {
    if (!has_parameter_classification(requirements)) {
        return result.classifications.empty();
    }

    return result.classifications.size() == molecules.size();
}

[[nodiscard]] auto make_invalid_classification_count_issue(
    const Method& method, const features::PreparedMoleculeCollection& molecules,
    const CollectionParameterPrerequisiteResult& result) -> PrerequisiteIssue {
    return make_issue(PrerequisiteIssueKind::parameter_classification_failed,
                      "method '" + std::string{method.id()} + "' produced " +
                          std::to_string(result.classifications.size()) +
                          " parameter classifications for " + std::to_string(molecules.size()) +
                          " molecules");
}

} // namespace

auto find_applicable_methods(const features::PreparedMoleculeCollection& molecules,
                             std::span<const Method* const> methods,
                             std::span<const parameters::ParameterSet> parameter_sets)
    -> ApplicabilityResult {
    ApplicabilityResult result;

    for (std::size_t method_index = 0; method_index < methods.size(); ++method_index) {
        const auto* method = methods[method_index];

        if (method == nullptr) {
            throw std::invalid_argument{"applicability check received a null method"};
        }

        MethodOptions method_options;

        try {
            method_options = make_default_options(method->option_schema());
        } catch (const std::exception& error) {
            add_rejected(result, method_index, std::nullopt,
                         {make_issue(PrerequisiteIssueKind::invalid_options, error.what())});

            continue;
        }

        const auto method_result = check_method_prerequisites(*method, molecules, method_options);

        if (!method_result) {
            add_rejected(result, method_index, std::nullopt, copy_issues(method_result.issues()));

            continue;
        }

        const auto requirements = method->requirements();

        if (!requirements.requires_parameters()) {
            result.applicable.push_back(ApplicableMethod{.method = method,
                                                         .parameter_set = nullptr,
                                                         .method_options = method_options,
                                                         .classifications = {}});

            continue;
        }

        if (parameter_sets.empty()) {
            add_rejected(
                result, method_index, std::nullopt,
                {make_issue(PrerequisiteIssueKind::missing_parameters,
                            "method '" + std::string{method->id()} +
                                "' requires parameters, but no parameter sets were provided")});

            continue;
        }

        for (std::size_t parameter_set_index = 0; parameter_set_index < parameter_sets.size();
             ++parameter_set_index) {
            const auto& parameter_set = parameter_sets[parameter_set_index];

            auto parameter_result =
                check_parameter_prerequisites(*method, molecules, parameter_set);

            if (!parameter_result) {
                add_rejected(result, method_index, parameter_set_index,
                             std::move(parameter_result.issues));

                continue;
            }

            if (!valid_classification_count(requirements, molecules, parameter_result)) {
                add_rejected(result, method_index, parameter_set_index,
                             {make_invalid_classification_count_issue(*method, molecules,
                                                                      parameter_result)});

                continue;
            }

            result.applicable.push_back(
                ApplicableMethod{.method = method,
                                 .parameter_set = &parameter_set,
                                 .method_options = method_options,
                                 .classifications = std::move(parameter_result.classifications)});
        }
    }

    return result;
}

} // namespace chargefw::methods