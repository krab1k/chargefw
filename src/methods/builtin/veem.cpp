#include "methods/builtin/veem.h"

#include <chargefw/core/periodic_table.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

[[nodiscard]] auto valence_electron_count_or_throw(const core::Element& element) -> int {
    const auto count = element.valence_electron_count();

    if (!count.has_value()) {
        throw std::logic_error{
            "VEEM calculation called for element '" + std::string{element.symbol} +
            "' without available valence electron count"
        };
    }

    return *count;
}

} // namespace

auto VEEMMethod::add_method_specific_prerequisite_issues(
    const MethodPrerequisiteInput& input,
    PrerequisiteResult& result
) const -> void {
    const auto& table = core::periodic_table();
    const auto& molecule = input.prepared_molecule.molecule();

    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto& atom = molecule.atom(atom_index);

        if (!table.contains(atom.atomic_number())) {
            result.add(PrerequisiteIssue{
                .kind = PrerequisiteIssueKind::unsupported_molecule,
                .message = "atom " + std::to_string(atom_index) +
                           " has atomic number outside the bundled periodic table",
                .atom_index = atom_index
            });
            continue;
        }

        const auto& element = table.element(atom.atomic_number());

        if (!element.valence_electron_count().has_value()) {
            result.add(PrerequisiteIssue{
                .kind = PrerequisiteIssueKind::unsupported_molecule,
                .message = "atom " + std::to_string(atom_index) + " (" +
                           std::string{element.symbol} +
                           ") has no VEEM valence electron count",
                .atom_index = atom_index
            });
        }
    }
}

auto VEEMMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& table = core::periodic_table();
    const auto& molecule = input.prepared_molecule.molecule();

    std::vector<double> charges(molecule.atom_count(), 0.0);

    if (molecule.atom_count() == 0) {
        return charges::AtomicCharges{std::move(charges)};
    }

    auto numerator = 0.0;
    auto denominator = 0.0;

    for (const auto& atom : molecule.atoms()) {
        const auto& element = table.element(atom.atomic_number());
        const auto valence = static_cast<double>(valence_electron_count_or_throw(element));

        numerator += element.electronegativity * valence;
        denominator += valence;
    }

    if (denominator == 0.0) {
        throw std::logic_error{"VEEM calculation has zero valence-electron denominator"};
    }

    const auto equalized_electronegativity = numerator / denominator;

    if (equalized_electronegativity == 0.0) {
        throw std::logic_error{"VEEM calculation has zero equalized electronegativity"};
    }

    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto& atom = molecule.atom(atom_index);
        const auto& element = table.element(atom.atomic_number());
        const auto valence = static_cast<double>(valence_electron_count_or_throw(element));

        charges[atom_index] =
            valence * (equalized_electronegativity - element.electronegativity) /
            equalized_electronegativity;
    }

    return charges::AtomicCharges{std::move(charges)};
}

} // namespace chargefw::methods::builtin