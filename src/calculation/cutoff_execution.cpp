#include "calculation/cutoff_execution.h"
#include "calculation/parallel_for.h"
#include "calculation/reduced_execution.h"
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
#include <vector>

namespace chargefw::calculation {
namespace {

[[nodiscard]] auto
calculate_target(const methods::ApplicableMethod& selected,
                 const features::PreparedMolecule& source,
                 const parameters::ParameterClassification* source_classification,
                 const std::size_t conformer_index, const double radius,
                 const ChargeCorrectionPolicy charge_correction, const std::size_t max_threads)
    -> charges::AtomicCharges {
    const auto& source_molecule = source.molecule();
    auto values = std::vector<double>(source_molecule.atom_count());
    const auto requirements = selected.method->requirements();
    const features::ConformerFeatures source_geometry{source_molecule, conformer_index};
    const features::SpatialFragmentBuilder fragment_builder{source, source_geometry};

    ::chargefw::calculation::detail::parallel_for_indexed(
        source_molecule.atom_count(), max_threads, [&](const std::size_t center_source_atom_index) {
            try {
                const auto fragment = fragment_builder.build(center_source_atom_index, radius);
                const auto fragment_charges = detail::calculate_fragment_charges(
                    selected, source_molecule, source_classification, fragment);
                values[center_source_atom_index] =
                    fragment_charges[fragment.center_local_atom_index()];
            } catch (const std::exception& error) {
                throw std::runtime_error{
                    "cutoff calculation failed for method '" + std::string{selected.method->id()} +
                    "', molecule '" + std::string{source_molecule.name()} + "', conformer " +
                    std::to_string(conformer_index) + ", center atom " +
                    std::to_string(center_source_atom_index) + ": " + error.what()};
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
                              const ExecutionPolicy& policy, const std::size_t max_threads)
    -> charges::ChargeSet {
    methods::detail::validate_selected_candidate(selected, molecules);
    methods::detail::validate_coordinate_targets(selected, molecules);
    detail::validate_reduced_request(selected, policy, ExecutionMode::cutoff);
    const auto radius = policy.radius();
    if (!radius.has_value()) {
        throw std::logic_error{"validated cutoff execution policy has no radius"};
    }

    struct Target {
        std::size_t molecule_index;
        std::size_t conformer_index;
    };

    std::vector<Target> targets;
    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        const auto& molecule = molecules[molecule_index];
        for (std::size_t conformer_index = 0;
             conformer_index < molecule.molecule().conformer_count(); ++conformer_index) {
            targets.push_back(
                Target{.molecule_index = molecule_index, .conformer_index = conformer_index});
        }
    }

    std::vector<std::optional<charges::ChargeAssignment>> calculated(targets.size());
    const auto target_threads = targets.size() > 1 ? max_threads : 1;
    ::chargefw::calculation::detail::parallel_for_indexed(
        targets.size(), target_threads, [&](const std::size_t target_index) {
            const auto& target = targets[target_index];
            const auto& molecule = molecules[target.molecule_index];
            const auto* classification =
                selected.uses_parameters()
                    ? std::addressof(selected.classifications[target.molecule_index])
                    : nullptr;
            calculated[target_index] = charges::ChargeAssignment{
                .target = charges::ChargeTarget{.molecule_index = target.molecule_index,
                                                .conformer_index = target.conformer_index},
                .charges = calculate_target(
                    selected, molecule, classification, target.conformer_index, *radius,
                    policy.charge_correction(), target_threads == 1 ? max_threads : 1)};
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
