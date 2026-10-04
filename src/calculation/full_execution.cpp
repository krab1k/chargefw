#include "calculation/full_execution.h"

#include "calculation/fixed_charge_partition.h"
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
        throw std::runtime_error{"method '" + std::string{method.id()} + "' produced " +
                                 std::to_string(atomic_charges.size()) + " charges for molecule '" +
                                 std::string{molecule.molecule().name()} + "' with " +
                                 std::to_string(expected_size) + " atoms"};
    }
}

[[nodiscard]] auto calculate_method(const methods::Method& method,
                                    const methods::CalculationInput& input)
    -> charges::AtomicCharges {
    try {
        return method.calculate(input);
    } catch (const std::invalid_argument& error) {
        throw std::runtime_error{error.what()};
    }
}

[[nodiscard]] auto calculate_target(const methods::ApplicableMethod& selected,
                                    const features::PreparedMolecule& molecule,
                                    const parameters::ParameterClassification* classification,
                                    const std::optional<std::size_t> conformer_index,
                                    const detail::FixedChargePartitionTarget* partition_target)
    -> charges::AtomicCharges {
    auto geometry = make_geometry_if_required(*selected.method, molecule, conformer_index);
    auto fixed_sources = std::vector<methods::FixedPointSource>{};
    auto target_charge = core::total_formal_charge(molecule.molecule());
    if (partition_target != nullptr) {
        target_charge = partition_target->active_charge;
        if (!partition_target->sources.empty()) {
            if (!conformer_index.has_value()) {
                throw std::invalid_argument{"fixed-charge embedding requires a conformer index"};
            }
            const auto& positions = partition_target->source_positions.at(*conformer_index);
            fixed_sources.reserve(partition_target->sources.size());
            for (std::size_t index = 0; index < partition_target->sources.size(); ++index) {
                fixed_sources.push_back(
                    {positions.at(index), partition_target->sources.at(index).charge});
            }
        }
    }

    if (classification == nullptr) {
        const methods::CalculationInput input{
            molecule,      selected.method_options,
            target_charge, geometry ? std::addressof(*geometry) : nullptr,
            nullptr,       fixed_sources};
        auto atomic_charges = calculate_method(*selected.method, input);
        validate_calculated_charges(*selected.method, molecule, atomic_charges);
        return atomic_charges;
    }

    parameters::validate_parameter_classification(molecule.molecule(), *selected.parameter_set,
                                                  *classification);
    const parameters::ParameterView parameter_view{*selected.parameter_set, *classification};
    const methods::CalculationInput input{
        molecule,        selected.method_options,
        target_charge,   geometry ? std::addressof(*geometry) : nullptr,
        &parameter_view, fixed_sources};
    auto atomic_charges = calculate_method(*selected.method, input);
    validate_calculated_charges(*selected.method, molecule, atomic_charges);
    return atomic_charges;
}

} // namespace

auto calculate_full_charges(const methods::ApplicableMethod& selected,
                            const features::PreparedMoleculeCollection& molecules,
                            const std::size_t max_threads, const CalculationObserver& observer,
                            const detail::FixedChargePartition* fixed_charge_partition)
    -> charges::ChargeSet {
    methods::detail::validate_selected_candidate(selected, molecules,
                                                 fixed_charge_partition != nullptr);
    if (fixed_charge_partition != nullptr) {
        if (fixed_charge_partition->active_molecules.size() != molecules.size() ||
            fixed_charge_partition->targets.size() != molecules.size()) {
            throw std::invalid_argument{
                "prepared molecule collection does not match fixed-charge partition"};
        }
        for (std::size_t index = 0; index < molecules.size(); ++index) {
            if (std::addressof(molecules[index].molecule()) !=
                std::addressof(fixed_charge_partition->active_molecules[index])) {
                throw std::invalid_argument{
                    "prepared molecule collection does not match fixed-charge partition"};
            }
        }
    }
    methods::detail::validate_coordinate_targets(selected, molecules);

    return detail::execute_calculation_targets(
        selected, molecules, ExecutionMode::full, selected.method->requirements().coordinates,
        detail::ParallelizationLevel::targets, max_threads, observer,
        [&](const features::PreparedMolecule& molecule,
            const parameters::ParameterClassification* classification,
            const std::optional<std::size_t> conformer_index, const std::size_t,
            const detail::ProgressContext& context) {
            const auto* partition_target =
                fixed_charge_partition == nullptr
                    ? nullptr
                    : std::addressof(fixed_charge_partition->targets.at(context.molecule_index));
            return calculate_target(selected, molecule, classification, conformer_index,
                                    partition_target);
        });
}

} // namespace chargefw::calculation
