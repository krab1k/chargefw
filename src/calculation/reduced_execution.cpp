#include "calculation/reduced_execution.h"

#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/models/parameter_view.h>

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
    if (charges.size() != fragment.molecule().atom_count()) {
        throw std::runtime_error{"method '" + std::string{selected.method->id()} + "' produced " +
                                 std::to_string(charges.size()) + " fragment charges for " +
                                 std::to_string(fragment.molecule().atom_count()) + " atoms"};
    }
    return charges;
}

} // namespace chargefw::calculation::detail
