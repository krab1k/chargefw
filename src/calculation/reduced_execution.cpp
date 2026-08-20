#include "calculation/reduced_execution.h"

#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::calculation::detail {

auto validate_reduced_request(const methods::ApplicableMethod& selected,
                              const ExecutionPolicy& policy, const ExecutionMode mode) -> void {
    if (policy.mode() != mode) {
        throw std::invalid_argument{"reduced executor received an unexpected execution policy"};
    }
    if (selected.method == nullptr) {
        throw std::invalid_argument{"selected applicable method has no method"};
    }

    const auto requirements = selected.method->requirements();
    const auto supports_mode = mode == ExecutionMode::cutoff
                                   ? requirements.resources.supports_cutoff
                                   : requirements.resources.supports_cover;
    if (!requirements.coordinates || !supports_mode) {
        throw std::invalid_argument{"method '" + std::string{selected.method->id()} +
                                    "' does not support " + std::string{to_string(mode)} +
                                    " execution"};
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

auto fragment_target_charge(const methods::FragmentTargetChargePolicy policy,
                            const core::Molecule& source, const features::SpatialFragment& fragment)
    -> double {
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

auto final_target_charge(const methods::FragmentTargetChargePolicy policy,
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

auto parameter_set_id_for(const methods::ApplicableMethod& selected) -> std::optional<std::string> {
    if (!selected.uses_parameters()) {
        return std::nullopt;
    }
    return std::string{selected.parameter_set->id()};
}

auto calculate_fragment_charges(const methods::ApplicableMethod& selected,
                                const core::Molecule& source,
                                const parameters::ParameterClassification* source_classification,
                                const features::SpatialFragment& fragment)
    -> charges::AtomicCharges {
    const features::PreparedMolecule prepared_fragment{fragment.molecule()};
    const features::ConformerFeatures geometry{fragment.molecule(), 0};
    const auto target_charge = fragment_target_charge(
        selected.method->requirements().resources.fragment_target_charge_policy, source, fragment);

    auto projected_classification = parameters::ParameterClassification{};
    if (source_classification != nullptr) {
        projected_classification =
            features::project_classification(*source_classification, fragment);
        parameters::validate_parameter_classification(fragment.molecule(), *selected.parameter_set,
                                                      projected_classification);
    }
    const auto parameter_view =
        source_classification == nullptr
            ? std::optional<parameters::ParameterView>{}
            : std::optional<parameters::ParameterView>{std::in_place, *selected.parameter_set,
                                                       projected_classification};

    const methods::CalculationInput input{
        prepared_fragment, selected.method_options, target_charge, &geometry,
        parameter_view ? std::addressof(*parameter_view) : nullptr};
    const auto charges = selected.method->calculate(input);
    validate_fragment_charges(*selected.method, fragment, charges);
    return charges;
}

} // namespace chargefw::calculation::detail
