#include "calculation/cover_execution.h"

#include "calculation/reduced_execution.h"
#include "calculation/target_execution.h"
#include "methods/applicable_method_execution.h"

#include <chargefw/features/conformer_features.h>
#include <chargefw/features/spatial_fragment.h>
#include <chargefw/methods/method.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::calculation {
namespace {

inline constexpr double cover_retained_radius = 3.0;

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

    struct WorkItem {
        std::size_t pivot_source_atom_index;
        std::size_t owned_count = 0;
    };

    const auto unassigned = std::numeric_limits<std::size_t>::max();
    auto owners = std::vector<std::size_t>(source_molecule.atom_count(), unassigned);
    auto work_items = std::vector<WorkItem>{};
    for (std::size_t pivot_source_atom_index = 0;
         pivot_source_atom_index < source_molecule.atom_count(); ++pivot_source_atom_index) {
        if (owners[pivot_source_atom_index] != unassigned) {
            continue;
        }

        const auto work_item_index = work_items.size();
        work_items.push_back(WorkItem{.pivot_source_atom_index = pivot_source_atom_index});
        for (std::size_t source_atom_index = 0; source_atom_index < source_molecule.atom_count();
             ++source_atom_index) {
            if (owners[source_atom_index] == unassigned &&
                source_geometry.distance(pivot_source_atom_index, source_atom_index) <=
                    cover_retained_radius) {
                owners[source_atom_index] = work_item_index;
                ++work_items.back().owned_count;
            }
        }

        if (owners[pivot_source_atom_index] != work_item_index) {
            throw std::logic_error{"cover pivot did not retain its source atom"};
        }
    }

    if (std::ranges::find(owners, unassigned) != owners.end()) {
        throw std::logic_error{"cover calculation left a source atom without an owner"};
    }

    ::chargefw::calculation::detail::progress_for_indexed(
        work_items.size(), max_threads, progress_ctx, [&](const std::size_t work_item_index) {
            const auto pivot_source_atom_index =
                work_items[work_item_index].pivot_source_atom_index;

            try {
                const auto fragment = fragment_builder.build(pivot_source_atom_index, radius);
                const auto fragment_charges = detail::calculate_fragment_charges(
                    selected, source_molecule, source_classification, fragment);
                const auto local_to_source = fragment.local_to_source_atom_indices();
                auto assigned_count = std::size_t{0};
                for (std::size_t local_atom_index = 0; local_atom_index < local_to_source.size();
                     ++local_atom_index) {
                    const auto source_atom_index = local_to_source[local_atom_index];
                    if (owners[source_atom_index] == work_item_index) {
                        values[source_atom_index] = fragment_charges[local_atom_index];
                        ++assigned_count;
                    }
                }
                if (assigned_count != work_items[work_item_index].owned_count) {
                    throw std::logic_error{
                        "cover fragment did not include every owned source atom"};
                }
            } catch (const std::exception& error) {
                throw std::runtime_error{
                    "cover calculation failed for method '" + std::string{selected.method->id()} +
                    "', molecule '" + std::string{source_molecule.name()} + "', conformer " +
                    std::to_string(conformer_index) + ", pivot atom " +
                    std::to_string(pivot_source_atom_index) + ": " + error.what()};
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

auto calculate_cover_charges(const methods::ApplicableMethod& selected,
                             const features::PreparedMoleculeCollection& molecules,
                             const ExecutionPolicy& policy, const std::size_t max_threads,
                             const CalculationObserver& observer) -> charges::ChargeSet {
    methods::detail::validate_selected_candidate(selected, molecules);
    methods::detail::validate_coordinate_targets(selected, molecules);
    detail::validate_reduced_request(selected, policy, ExecutionMode::cover);
    const auto radius = policy.radius();
    if (!radius.has_value()) {
        throw std::logic_error{"validated cover execution policy has no radius"};
    }
    return detail::execute_calculation_targets(
        selected, molecules, ExecutionMode::cover, true, detail::ParallelizationLevel::fragments,
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
