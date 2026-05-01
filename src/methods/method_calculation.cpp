#include <chargefw/methods/method_calculation.h>

#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/parameter_view.h>

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::methods {
namespace {

auto validate_selected_candidate(const ApplicableMethod& selected,
                                 const features::PreparedMoleculeCollection& molecules) -> void {
    if (selected.method == nullptr) {
        throw std::invalid_argument{"selected applicable method has no method"};
    }

    if (selected.method->requires_parameters() && !selected.uses_parameters()) {
        throw std::invalid_argument{"selected method '" + std::string{selected.method->id()} +
                                    "' requires parameters, but no parameter set is attached"};
    }

    if (!selected.method->requires_parameters() && selected.uses_parameters()) {
        throw std::invalid_argument{
            "selected method '" + std::string{selected.method->id()} +
            "' does not require parameters, but a parameter set is attached"};
    }

    if (!selected.uses_parameters()) {
        if (!selected.classifications.empty()) {
            throw std::invalid_argument{"selected method '" + std::string{selected.method->id()} +
                                        "' has classifications, but no parameter set"};
        }

        return;
    }

    if (selected.classifications.size() != molecules.molecule_count()) {
        throw std::invalid_argument{"selected method '" + std::string{selected.method->id()} +
                                    "' has " + std::to_string(selected.classifications.size()) +
                                    " classifications for " +
                                    std::to_string(molecules.molecule_count()) + " molecules"};
    }
}

[[nodiscard]] auto parameter_set_id_for(const ApplicableMethod& selected)
    -> std::optional<std::string> {
    if (!selected.uses_parameters()) {
        return std::nullopt;
    }

    return std::string{selected.parameter_set->id()};
}

[[nodiscard]] auto calculate_without_parameters(const ApplicableMethod& selected,
                                                const features::PreparedMolecule& molecule)
    -> charges::AtomicCharges {
    const CalculationInput input{molecule, selected.method_options};

    return selected.method->calculate(input);
}

[[nodiscard]] auto calculate_with_parameters(
    const ApplicableMethod& selected, const features::PreparedMolecule& molecule,
    const parameters::ParameterClassification& classification) -> charges::AtomicCharges {
    const parameters::ParameterView parameter_view{*selected.parameter_set, classification};

    const CalculationInput input{molecule, selected.method_options, nullptr, &parameter_view};

    return selected.method->calculate(input);
}

} // namespace

auto calculate_charges(const ApplicableMethod& selected,
                       const features::PreparedMoleculeCollection& molecules)
    -> charges::ChargeSet {
    validate_selected_candidate(selected, molecules);

    std::vector<charges::ChargeAssignment> assignments;
    assignments.reserve(molecules.molecule_count());

    for (std::size_t molecule_index = 0; molecule_index < molecules.molecule_count();
         ++molecule_index) {
        auto atomic_charges =
            selected.uses_parameters()
                ? calculate_with_parameters(selected, molecules[molecule_index],
                                            selected.classifications[molecule_index])
                : calculate_without_parameters(selected, molecules[molecule_index]);

        assignments.push_back(charges::ChargeAssignment{
            .target = charges::ChargeTarget{.molecule_index = molecule_index,
                                            .conformer_index = std::nullopt},
            .charges = std::move(atomic_charges)});
    }

    return charges::ChargeSet{std::string{selected.method->id()}, std::move(assignments),
                              parameter_set_id_for(selected)};
}

} // namespace chargefw::methods