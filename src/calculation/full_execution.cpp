#include "calculation/full_execution.h"

#include "calculation/parallel_for.h"
#include "methods/applicable_method_execution.h"

#include <chargefw/core/molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <chrono>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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

    const auto geometry_dependent = selected.method->requirements().coordinates;
    struct Target {
        std::size_t molecule_index;
        std::optional<std::size_t> conformer_index;
    };

    std::vector<Target> targets;
    std::size_t target_count = 0;
    for (const auto& prepared_molecule : molecules.molecules()) {
        target_count += geometry_dependent ? prepared_molecule.molecule().conformer_count() : 1;
    }
    targets.reserve(target_count);

    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        const auto& molecule = molecules[molecule_index].molecule();
        const auto molecule_target_count = geometry_dependent ? molecule.conformer_count() : 1;
        for (std::size_t target_index = 0; target_index < molecule_target_count; ++target_index) {
            targets.push_back(Target{
                .molecule_index = molecule_index,
                .conformer_index =
                    geometry_dependent ? std::optional<std::size_t>{target_index} : std::nullopt});
        }
    }

    const auto computation_start = std::chrono::steady_clock::now();
    std::vector<std::optional<charges::ChargeAssignment>> calculated(targets.size());
    detail::parallel_for_indexed(targets.size(), max_threads, [&](const std::size_t target_index) {
        const auto& target = targets[target_index];
        const auto& molecule = molecules[target.molecule_index];
        const detail::ProgressContext target_ctx{
            .observer = observer,
            .mode = ExecutionMode::full,
            .method_id = selected.method->id(),
            .target_index = target_index,
            .target_count = targets.size(),
            .molecule_index = target.molecule_index,
            .conformer_index = target.conformer_index,
            .computation_start = computation_start,
        };

        detail::check_cancellation(observer);
        detail::emit_target_event(observer, CalculationPhase::target_started, target_ctx);
        detail::check_cancellation(observer);

        const auto* classification =
            selected.uses_parameters()
                ? std::addressof(selected.classifications[target.molecule_index])
                : nullptr;
        calculated[target_index] = charges::ChargeAssignment{
            .target = charges::ChargeTarget{.molecule_index = target.molecule_index,
                                            .conformer_index = target.conformer_index},
            .charges =
                calculate_target(selected, molecule, classification, target.conformer_index)};

        detail::emit_target_event(observer, CalculationPhase::target_finished, target_ctx);
    });

    std::vector<charges::ChargeAssignment> assignments;
    assignments.reserve(calculated.size());
    for (auto& assignment : calculated) {
        assignments.push_back(std::move(*assignment));
    }

    return charges::ChargeSet{std::string{selected.method->id()}, std::move(assignments),
                              methods::detail::parameter_set_id_for(selected)};
}

} // namespace chargefw::calculation
