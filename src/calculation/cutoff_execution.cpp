#include "calculation/cutoff_execution.h"
#include "calculation/reduced_execution.h"
#include "calculation/target_execution.h"
#include "methods/applicable_method_execution.h"

#include <chargefw/core/molecule.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/spatial_fragment.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/classification/parameter_classification.h>

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::calculation {
namespace {

[[nodiscard]] auto
calculate_target(const methods::ApplicableMethod& selected,
                 const features::PreparedMolecule& source,
                 const parameters::ParameterClassification* source_classification,
                 const std::size_t conformer_index, const double radius,
                 const ChargeCorrectionPolicy charge_correction, const std::size_t max_threads,
                 const detail::ProgressContext& progress_ctx) -> charges::AtomicCharges {
    const auto& source_molecule = source.molecule();
    auto values = std::vector<double>(source_molecule.atom_count());
    const auto requirements = selected.method->requirements();
    const features::ConformerFeatures source_geometry{source_molecule, conformer_index};
    const features::SpatialFragmentBuilder fragment_builder{source, source_geometry};

    ::chargefw::calculation::detail::progress_for_indexed(
        source_molecule.atom_count(), max_threads, progress_ctx,
        [&](const std::size_t center_source_atom_index) {
            try {
                const auto fragment = fragment_builder.build(center_source_atom_index, radius);
                const auto fragment_charges = detail::calculate_fragment_charges(
                    selected, source_molecule, source_classification, fragment);
                values[center_source_atom_index] =
                    fragment_charges[fragment.center_local_atom_index()];
            } catch (const std::exception& error) {
                throw std::runtime_error{"cutoff fragment around source atom " +
                                         std::to_string(center_source_atom_index + 1) +
                                         " failed: " + error.what()};
            }
        });

    detail::apply_charge_correction(
        values,
        detail::final_target_charge(requirements.resources.fragment_target_charge_policy,
                                    source_molecule),
        charge_correction);
    return charges::AtomicCharges{std::move(values)};
}

} // namespace

auto calculate_cutoff_charges(const methods::ApplicableMethod& selected,
                              const features::PreparedMoleculeCollection& molecules,
                              const ExecutionPolicy& policy, const std::size_t max_threads,
                              const CalculationObserver& observer) -> charges::ChargeSet {
    methods::detail::validate_selected_candidate(selected, molecules);
    methods::detail::validate_coordinate_targets(selected, molecules);
    detail::validate_reduced_request(selected, policy, ExecutionMode::cutoff);
    const auto radius = policy.radius();
    if (!radius.has_value()) {
        throw std::logic_error{"validated cutoff execution policy has no radius"};
    }

    return detail::execute_calculation_targets(
        selected, molecules, ExecutionMode::cutoff, true, detail::ParallelizationLevel::fragments,
        max_threads, observer,
        [&](const features::PreparedMolecule& molecule,
            const parameters::ParameterClassification* classification,
            const std::optional<std::size_t> conformer_index,
            const std::size_t fragment_max_threads, const detail::ProgressContext& target_ctx) {
            return calculate_target(selected, molecule, classification, *conformer_index, *radius,
                                    policy.charge_correction(), fragment_max_threads, target_ctx);
        });
}

} // namespace chargefw::calculation
