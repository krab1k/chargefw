#include "calculation/reduced_execution.h"
#include "features/topology_helpers.h"

#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <algorithm>
#include <memory>
#include <numeric>
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
    if (requirements.resources.reduced_charge_policy == methods::ReducedChargePolicy::unsupported) {
        throw std::invalid_argument{"method '" + std::string{selected.method->id()} +
                                    "' has no fragment target-charge policy"};
    }
}

auto prepare_reduced_charge_context(
    const methods::ApplicableMethod& selected, const features::PreparedMolecule& source,
    const parameters::ParameterClassification* source_classification) -> ReducedChargeContext {
    const auto& molecule = source.molecule();
    const auto policy = selected.method->requirements().resources.reduced_charge_policy;
    auto context = ReducedChargeContext{};
    context.reference_charges.resize(molecule.atom_count());

    switch (policy) {
    case methods::ReducedChargePolicy::uniform_target_global:
        if (!context.reference_charges.empty()) {
            const auto reference = core::total_formal_charge(molecule) /
                                   static_cast<double>(context.reference_charges.size());
            std::ranges::fill(context.reference_charges, reference);
        }
        if (!context.reference_charges.empty()) {
            context.conservation_groups.emplace_back(context.reference_charges.size());
            std::iota(context.conservation_groups.front().begin(),
                      context.conservation_groups.front().end(), std::size_t{0});
        }
        return context;
    case methods::ReducedChargePolicy::zero_components:
        break;
    case methods::ReducedChargePolicy::formal_charge_components:
        for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
            context.reference_charges[atom_index] = molecule.atom(atom_index).formal_charge();
        }
        break;
    case methods::ReducedChargePolicy::parameterized_charge_components: {
        if (selected.parameter_set == nullptr || source_classification == nullptr) {
            throw std::invalid_argument{
                "parameterized reduced charge policy requires classified parameters"};
        }
        const parameters::ParameterView parameters{*selected.parameter_set, *source_classification};
        const auto q0 = parameters.atom("q0");
        for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
            context.reference_charges[atom_index] = q0[atom_index];
        }
        if (!context.reference_charges.empty()) {
            const auto raw_total = std::accumulate(context.reference_charges.begin(),
                                                   context.reference_charges.end(), 0.0);
            const auto offset = (core::total_formal_charge(molecule) - raw_total) /
                                static_cast<double>(context.reference_charges.size());
            for (auto& reference : context.reference_charges) {
                reference += offset;
            }
        }
        break;
    }
    case methods::ReducedChargePolicy::unsupported:
        break;
    }

    if (policy == methods::ReducedChargePolicy::unsupported) {
        throw std::invalid_argument{"unsupported reduced charge policy"};
    }
    context.conservation_groups = features::connected_components(source.topology().adjacency());
    return context;
}

auto fragment_target_charge(const ReducedChargeContext& context,
                            const features::SpatialFragment& fragment) -> double {
    auto target = 0.0;
    for (const auto source_atom_index : fragment.local_to_source_atom_indices()) {
        target += context.reference_charges.at(source_atom_index);
    }
    return target;
}

auto enforce_conserved_charges(std::vector<double>& values, const ReducedChargeContext& context)
    -> void {
    if (values.size() != context.reference_charges.size()) {
        throw std::invalid_argument{"assembled charges and source references have different sizes"};
    }
    for (const auto& group : context.conservation_groups) {
        auto reference_total = 0.0;
        auto calculated_total = 0.0;
        for (const auto atom_index : group) {
            reference_total += context.reference_charges.at(atom_index);
            calculated_total += values.at(atom_index);
        }
        const auto delta = (reference_total - calculated_total) / static_cast<double>(group.size());
        for (const auto atom_index : group) {
            values[atom_index] += delta;
        }
    }
}

auto calculate_fragment_charges(const methods::ApplicableMethod& selected,
                                const parameters::ParameterClassification* source_classification,
                                const features::SpatialFragment& fragment,
                                const ReducedChargeContext& charge_context)
    -> charges::AtomicCharges {
    const features::PreparedMolecule prepared_fragment{fragment.molecule()};
    const features::ConformerFeatures geometry{fragment.molecule(), 0};
    const auto target_charge = fragment_target_charge(charge_context, fragment);

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
