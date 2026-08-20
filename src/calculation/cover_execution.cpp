#include "calculation/cover_execution.h"

#include "calculation/reduced_execution.h"

#include <chargefw/features/conformer_features.h>
#include <chargefw/features/spatial_fragment.h>
#include <chargefw/methods/method.h>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::calculation {
namespace {

inline constexpr double cover_retained_radius = 3.0;

[[nodiscard]] auto
calculate_target(const methods::ApplicableMethod& selected,
                 const features::PreparedMolecule& source,
                 const parameters::ParameterClassification* source_classification,
                 const std::size_t conformer_index, const double radius,
                 const ChargeCorrectionPolicy charge_correction) -> charges::AtomicCharges {
    const auto& source_molecule = source.molecule();
    auto values = std::vector<double>(source_molecule.atom_count());
    auto assigned = std::vector<bool>(source_molecule.atom_count(), false);
    const auto requirements = selected.method->requirements();
    const features::ConformerFeatures source_geometry{source_molecule, conformer_index};
    const features::SpatialFragmentBuilder fragment_builder{source, source_geometry};

    for (std::size_t pivot_source_atom_index = 0;
         pivot_source_atom_index < source_molecule.atom_count(); ++pivot_source_atom_index) {
        if (assigned[pivot_source_atom_index]) {
            continue;
        }

        try {
            const auto fragment = fragment_builder.build(pivot_source_atom_index, radius);
            const auto fragment_charges = detail::calculate_fragment_charges(
                selected, source_molecule, source_classification, fragment);
            const auto local_to_source = fragment.local_to_source_atom_indices();
            auto retained_count = std::size_t{0};
            for (std::size_t local_atom_index = 0; local_atom_index < local_to_source.size();
                 ++local_atom_index) {
                const auto source_atom_index = local_to_source[local_atom_index];
                if (assigned[source_atom_index] ||
                    source_geometry.distance(pivot_source_atom_index, source_atom_index) >
                        cover_retained_radius) {
                    continue;
                }

                values[source_atom_index] = fragment_charges[local_atom_index];
                assigned[source_atom_index] = true;
                ++retained_count;
            }

            if (retained_count == 0) {
                throw std::logic_error{"cover pivot did not retain its source atom"};
            }
        } catch (const std::exception& error) {
            throw std::runtime_error{"cover calculation failed for method '" +
                                     std::string{selected.method->id()} + "', molecule '" +
                                     std::string{source_molecule.name()} + "', conformer " +
                                     std::to_string(conformer_index) + ", pivot atom " +
                                     std::to_string(pivot_source_atom_index) + ": " + error.what()};
        }
    }

    if (std::ranges::find(assigned, false) != assigned.end()) {
        throw std::logic_error{"cover calculation left a source atom without an owner"};
    }

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
                             const ExecutionPolicy& policy) -> charges::ChargeSet {
    detail::validate_reduced_request(selected, policy, ExecutionMode::cover);
    const auto radius = policy.radius();
    if (!radius.has_value()) {
        throw std::logic_error{"validated cover execution policy has no radius"};
    }
    if (selected.uses_parameters() && selected.classifications.size() != molecules.size()) {
        throw std::invalid_argument{"selected method '" + std::string{selected.method->id()} +
                                    "' has an invalid number of classifications"};
    }

    auto assignments = std::vector<charges::ChargeAssignment>{};
    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        const auto& molecule = molecules[molecule_index];
        const auto* classification = selected.uses_parameters()
                                         ? std::addressof(selected.classifications[molecule_index])
                                         : nullptr;
        for (std::size_t conformer_index = 0;
             conformer_index < molecule.molecule().conformer_count(); ++conformer_index) {
            assignments.push_back(charges::ChargeAssignment{
                .target = charges::ChargeTarget{.molecule_index = molecule_index,
                                                .conformer_index = conformer_index},
                .charges = calculate_target(selected, molecule, classification, conformer_index,
                                            *radius, policy.charge_correction())});
        }
    }

    return charges::ChargeSet{std::string{selected.method->id()}, std::move(assignments),
                              detail::parameter_set_id_for(selected)};
}

} // namespace chargefw::calculation
