#include <chargefw/methods/method.h>

#include <cmath>
#include <exception>
#include <map>
#include <string>

namespace chargefw::methods {
namespace {

auto add_invalid_options_issue(PrerequisiteResult& result, const std::exception& error) -> void {
    result.add(
        PrerequisiteIssue{.kind = PrerequisiteIssueKind::invalid_options, .message = error.what()});
}

using PositionKey = core::Position;

struct PositionKeyLess {
    [[nodiscard]] auto operator()(const core::Position& first,
                                  const core::Position& second) const noexcept -> bool {
        const auto normalize = [](const double value) noexcept -> double {
            return value == 0.0 ? 0.0 : value;
        };

        const auto first_x = normalize(first.x);
        const auto first_y = normalize(first.y);
        const auto first_z = normalize(first.z);
        const auto second_x = normalize(second.x);
        const auto second_y = normalize(second.y);
        const auto second_z = normalize(second.z);

        if (first_x != second_x) {
            return first_x < second_x;
        }
        if (first_y != second_y) {
            return first_y < second_y;
        }
        return first_z < second_z;
    }
};

auto add_coordinate_prerequisite_issues(const Method& method, const core::Molecule& molecule,
                                        PrerequisiteResult& result) -> void {
    for (std::size_t conformer_index = 0; conformer_index < molecule.conformer_count();
         ++conformer_index) {
        const auto positions = molecule.conformer(conformer_index).positions();
        std::map<PositionKey, std::size_t, PositionKeyLess> first_atom_by_position;

        for (std::size_t atom_index = 0; atom_index < positions.size(); ++atom_index) {
            const auto& position = positions[atom_index];
            if (!std::isfinite(position.x) || !std::isfinite(position.y) ||
                !std::isfinite(position.z)) {
                result.add(
                    PrerequisiteIssue{.kind = PrerequisiteIssueKind::invalid_geometry,
                                      .message = "method '" + std::string{method.id()} +
                                                 "' requires finite coordinates; conformer " +
                                                 std::to_string(conformer_index) + ", atom " +
                                                 std::to_string(atom_index) + " is non-finite",
                                      .atom_index = atom_index});
                break;
            }

            const auto [found, inserted] = first_atom_by_position.emplace(position, atom_index);
            if (!inserted) {
                result.add(PrerequisiteIssue{
                    .kind = PrerequisiteIssueKind::invalid_geometry,
                    .message = "method '" + std::string{method.id()} +
                               "' requires distinct atom coordinates; conformer " +
                               std::to_string(conformer_index) + ", atoms " +
                               std::to_string(found->second) + " and " +
                               std::to_string(atom_index) + " share coordinates",
                    .atom_index = atom_index});
                break;
            }
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
