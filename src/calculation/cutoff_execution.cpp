#include "calculation/cutoff_execution.h"

#include <chargefw/core/molecule.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/spatial_fragment.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::calculation {
namespace {

auto validate_cutoff_request(const methods::ApplicableMethod& selected,
                             const ExecutionPolicy& policy) -> void {
    if (policy.mode() != ExecutionMode::cutoff) {
        throw std::invalid_argument{"cutoff executor requires a cutoff execution policy"};
    }
    if (selected.method == nullptr) {
        throw std::invalid_argument{"selected applicable method has no method"};
    }

    const auto requirements = selected.method->requirements();
    if (!requirements.coordinates || !requirements.resources.supports_cutoff) {
        throw std::invalid_argument{"method '" + std::string{selected.method->id()} +
                                    "' does not support cutoff execution"};
    }
    if (requirements.resources.fragment_target_charge_policy ==
        methods::FragmentTargetChargePolicy::unsupported) {
        throw std::invalid_argument{"method '" + std::string{selected.method->id()} +
                                    "' has no fragment target-charge policy"};
    }
    if (selected.method->requires_parameters() && !selected.uses_parameters()) {
        throw std::invalid_argument{"selected method '" + std::string{selected.method->id()} +
                                    "' requires parameters, but no parameter set is attached"};
    }
}

[[nodiscard]] auto fragment_target_charge(const methods::FragmentTargetChargePolicy policy,
                                          const core::Molecule& source,
                                          const features::SpatialFragment& fragment) -> double {
    switch (policy) {
    case methods::FragmentTargetChargePolicy::zero:
        return 0.0;
    case methods::FragmentTargetChargePolicy::proportional_to_atom_count:
        return static_cast<double>(fragment.molecule().atom_count()) /
               static_cast<double>(source.atom_count()) * core::total_formal_charge(source);
    case methods::FragmentTargetChargePolicy::unsupported:
        break;
    }

    throw std::invalid_argument{"unsupported fragment target-charge policy"};
}

[[nodiscard]] auto final_target_charge(const methods::FragmentTargetChargePolicy policy,
                                       const core::Molecule& source) -> double {
    switch (policy) {
    case methods::FragmentTargetChargePolicy::zero:
        return 0.0;
    case methods::FragmentTargetChargePolicy::proportional_to_atom_count:
        return core::total_formal_charge(source);
    case methods::FragmentTargetChargePolicy::unsupported:
        break;
    }

    throw std::invalid_argument{"unsupported fragment target-charge policy"};
}

auto validate_fragment_charges(const methods::Method& method,
                               const features::SpatialFragment& fragment,
                               const charges::AtomicCharges& charges) -> void {
    if (charges.size() != fragment.molecule().atom_count()) {
        throw std::runtime_error{"method '" + std::string{method.id()} + "' produced " +
                                 std::to_string(charges.size()) + " fragment charges for " +
                                 std::to_string(fragment.molecule().atom_count()) + " atoms"};
    }

    for (const auto value : charges.values()) {
        if (!std::isfinite(value)) {
            throw std::runtime_error{"method '" + std::string{method.id()} +
                                     "' produced a non-finite fragment charge"};
        }
    }
}

auto apply_charge_correction(std::vector<double>& values, const double target_charge,
                             const ChargeCorrectionPolicy policy) -> void {
    if (policy == ChargeCorrectionPolicy::none || values.empty()) {
        return;
    }

    if (policy != ChargeCorrectionPolicy::uniform) {
        throw std::invalid_argument{"unsupported charge correction policy"};
    }

    double total = 0.0;
    for (const auto value : values) {
        total += value;
    }
    const auto delta = (target_charge - total) / static_cast<double>(values.size());
    for (auto& value : values) {
        value += delta;
    }
}

[[nodiscard]] auto parameter_set_id_for(const methods::ApplicableMethod& selected)
    -> std::optional<std::string> {
    if (!selected.uses_parameters()) {
        return std::nullopt;
    }
    return std::string{selected.parameter_set->id()};
}

[[nodiscard]] auto
calculate_target(const methods::ApplicableMethod& selected,
                 const features::PreparedMolecule& source,
                 const parameters::ParameterClassification* source_classification,
                 const std::size_t conformer_index, const double radius,
                 const ChargeCorrectionPolicy charge_correction) -> charges::AtomicCharges {
    const auto& source_molecule = source.molecule();
    auto values = std::vector<double>(source_molecule.atom_count());
    const auto requirements = selected.method->requirements();
    const features::ConformerFeatures source_geometry{source_molecule, conformer_index};
    const features::SpatialFragmentBuilder fragment_builder{source, source_geometry};

    for (std::size_t center_source_atom_index = 0;
         center_source_atom_index < source_molecule.atom_count(); ++center_source_atom_index) {
        try {
            const auto fragment = fragment_builder.build(center_source_atom_index, radius);
            const features::PreparedMolecule prepared_fragment{fragment.molecule()};
            const features::ConformerFeatures geometry{fragment.molecule(), 0};
            const auto target_charge = fragment_target_charge(
                requirements.resources.fragment_target_charge_policy, source_molecule, fragment);

            auto projected_classification = parameters::ParameterClassification{};
            if (source_classification != nullptr) {
                projected_classification =
                    features::project_classification(*source_classification, fragment);
                parameters::validate_parameter_classification(
                    fragment.molecule(), *selected.parameter_set, projected_classification);
            }
            const auto parameter_view =
                source_classification == nullptr
                    ? std::optional<parameters::ParameterView>{}
                    : std::optional<parameters::ParameterView>{
                          std::in_place, *selected.parameter_set, projected_classification};

            const methods::CalculationInput input{
                prepared_fragment, selected.method_options, target_charge, &geometry,
                parameter_view ? std::addressof(*parameter_view) : nullptr};
            const auto fragment_charges = selected.method->calculate(input);
            validate_fragment_charges(*selected.method, fragment, fragment_charges);
            values[center_source_atom_index] = fragment_charges[fragment.center_local_atom_index()];
        } catch (const std::exception& error) {
            throw std::runtime_error{
                "cutoff calculation failed for method '" + std::string{selected.method->id()} +
                "', molecule '" + std::string{source_molecule.name()} + "', conformer " +
                std::to_string(conformer_index) + ", center atom " +
                std::to_string(center_source_atom_index) + ": " + error.what()};
        }
    }

    apply_charge_correction(
        values,
        final_target_charge(requirements.resources.fragment_target_charge_policy, source_molecule),
        charge_correction);
    return charges::AtomicCharges{std::move(values)};
}

} // namespace

auto calculate_cutoff_charges(const methods::ApplicableMethod& selected,
                              const features::PreparedMoleculeCollection& molecules,
                              const ExecutionPolicy& policy) -> charges::ChargeSet {
    validate_cutoff_request(selected, policy);
    const auto radius = policy.radius();
    if (!radius.has_value()) {
        throw std::logic_error{"validated cutoff execution policy has no radius"};
    }

    if (selected.uses_parameters() && selected.classifications.size() != molecules.size()) {
        throw std::invalid_argument{"selected method '" + std::string{selected.method->id()} +
                                    "' has an invalid number of classifications"};
    }

    std::vector<charges::ChargeAssignment> assignments;
    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        const auto& molecule = molecules[molecule_index];
        for (std::size_t conformer_index = 0;
             conformer_index < molecule.molecule().conformer_count(); ++conformer_index) {
            const auto* classification =
                selected.uses_parameters()
                    ? std::addressof(selected.classifications[molecule_index])
                    : nullptr;
            assignments.push_back(charges::ChargeAssignment{
                .target = charges::ChargeTarget{.molecule_index = molecule_index,
                                                .conformer_index = conformer_index},
                .charges = calculate_target(selected, molecule, classification, conformer_index,
                                            *radius, policy.charge_correction())});
        }
    }

    return charges::ChargeSet{std::string{selected.method->id()}, std::move(assignments),
                              parameter_set_id_for(selected)};
}

} // namespace chargefw::calculation
