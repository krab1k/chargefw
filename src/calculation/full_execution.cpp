#include "calculation/full_execution.h"

#include "calculation/target_execution.h"
#include "methods/applicable_method_execution.h"

#include <chargefw/core/molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::calculation {
namespace {

[[nodiscard]] auto make_geometry_if_required(const methods::Method& method,
                                             const features::PreparedMolecule& molecule,
                                             const std::optional<std::size_t> conformer_index)
    -> std::optional<features::ConformerFeatures> {
    if (!method.requirements().coordinates) {
        return std::nullopt;
    }

    if (!conformer_index.has_value()) {
        throw std::invalid_argument{"geometry-dependent method requires a conformer index"};
    }

    return features::ConformerFeatures{molecule.molecule(), *conformer_index};
}

auto validate_calculated_charges(const methods::Method& method,
                                 const features::PreparedMolecule& molecule,
                                 const charges::AtomicCharges& atomic_charges) -> void {
    const auto expected_size = molecule.molecule().atom_count();

    if (atomic_charges.size() != expected_size) {
        throw std::invalid_argument{"method '" + std::string{method.id()} + "' produced " +
                                    std::to_string(atomic_charges.size()) +
                                    " charges for molecule '" +
                                    std::string{molecule.molecule().name()} + "' with " +
                                    std::to_string(expected_size) + " atoms"};
    }
}

[[nodiscard]] auto calculate_target(const methods::ApplicableMethod& selected,
                                    const features::PreparedMolecule& molecule,
                                    const parameters::ParameterClassification* classification,
                                    const std::optional<std::size_t> conformer_index)
    -> charges::AtomicCharges {
    auto geometry = make_geometry_if_required(*selected.method, molecule, conformer_index);

    if (classification == nullptr) {
        const methods::CalculationInput input{molecule, selected.method_options,
                                              core::total_formal_charge(molecule.molecule()),
                                              geometry ? std::addressof(*geometry) : nullptr};
        auto atomic_charges = selected.method->calculate(input);
        validate_calculated_charges(*selected.method, molecule, atomic_charges);
        return atomic_charges;
    }

    parameters::validate_parameter_classification(molecule.molecule(), *selected.parameter_set,
                                                  *classification);
    const parameters::ParameterView parameter_view{*selected.parameter_set, *classification};
    const methods::CalculationInput input{
        molecule, selected.method_options, core::total_formal_charge(molecule.molecule()),
        geometry ? std::addressof(*geometry) : nullptr, &parameter_view};
    auto atomic_charges = selected.method->calculate(input);
    validate_calculated_charges(*selected.method, molecule, atomic_charges);
    return atomic_charges;
}

} // namespace

auto calculate_full_charges(const methods::ApplicableMethod& selected,
                            const features::PreparedMoleculeCollection& molecules,
                            const std::size_t max_threads, const CalculationObserver& observer)
    -> charges::ChargeSet {
    methods::detail::validate_selected_candidate(selected, molecules);
    methods::detail::validate_coordinate_targets(selected, molecules);

    return detail::execute_calculation_targets(
        selected, molecules, ExecutionMode::full, selected.method->requirements().coordinates,
        detail::ParallelizationLevel::targets, max_threads, observer,
        [&](const features::PreparedMolecule& molecule,
            const parameters::ParameterClassification* classification,
            const std::optional<std::size_t> conformer_index, const std::size_t,
            const detail::ProgressContext&) {
            return calculate_target(selected, molecule, classification, conformer_index);
        });
}

} // namespace chargefw::calculation
