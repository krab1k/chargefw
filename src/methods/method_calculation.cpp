#include <chargefw/methods/method_calculation.h>

#include "calculation/selected_candidate_validation.h"

#include <chargefw/core/molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_view.h>

#include "calculation/parallel_for.h"

#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chargefw::methods {
namespace {

[[nodiscard]] auto make_geometry_if_required(const Method& method,
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

auto validate_calculated_charges(const Method& method, const features::PreparedMolecule& molecule,
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

[[nodiscard]] auto parameter_set_id_for(const ApplicableMethod& selected)
    -> std::optional<std::string> {
    if (!selected.uses_parameters()) {
        return std::nullopt;
    }

    return std::string{selected.parameter_set->id()};
}

auto validate_coordinate_targets(const ApplicableMethod& selected,
                                 const features::PreparedMoleculeCollection& prepared_collection)
    -> void {
    if (!selected.method->requirements().coordinates) {
        return;
    }

    for (std::size_t molecule_index = 0; molecule_index < prepared_collection.size();
         ++molecule_index) {
        const auto& molecule = prepared_collection[molecule_index].molecule();

        if (molecule.conformer_count() == 0) {
            throw std::invalid_argument{"selected method '" + std::string{selected.method->id()} +
                                        "' requires coordinates, but molecule " +
                                        std::to_string(molecule_index) + " has no conformers"};
        }
    }
}

[[nodiscard]] auto calculate_without_parameters(const ApplicableMethod& selected,
                                                const features::PreparedMolecule& molecule,
                                                const std::optional<std::size_t> conformer_index)
    -> charges::AtomicCharges {
    auto geometry = make_geometry_if_required(*selected.method, molecule, conformer_index);

    const CalculationInput input{molecule, selected.method_options,
                                 core::total_formal_charge(molecule.molecule()),
                                 geometry ? std::addressof(*geometry) : nullptr};

    auto atomic_charges = selected.method->calculate(input);

    validate_calculated_charges(*selected.method, molecule, atomic_charges);

    return atomic_charges;
}

[[nodiscard]] auto calculate_with_parameters(
    const ApplicableMethod& selected, const features::PreparedMolecule& molecule,
    const parameters::ParameterClassification& classification,
    const std::optional<std::size_t> conformer_index) -> charges::AtomicCharges {
    parameters::validate_parameter_classification(molecule.molecule(), *selected.parameter_set,
                                                  classification);

    const parameters::ParameterView parameter_view{*selected.parameter_set, classification};

    auto geometry = make_geometry_if_required(*selected.method, molecule, conformer_index);

    const CalculationInput input{molecule, selected.method_options,
                                 core::total_formal_charge(molecule.molecule()),
                                 geometry ? std::addressof(*geometry) : nullptr, &parameter_view};

    auto atomic_charges = selected.method->calculate(input);

    validate_calculated_charges(*selected.method, molecule, atomic_charges);

    return atomic_charges;
}

} // namespace

auto calculate_charges(const ApplicableMethod& selected,
                       const features::PreparedMoleculeCollection& prepared_collection,
                       const std::size_t max_threads) -> charges::ChargeSet {
    calculation::detail::validate_selected_candidate(selected, prepared_collection);
    validate_coordinate_targets(selected, prepared_collection);

    const auto geometry_dependent = selected.method->requirements().coordinates;

    struct Target {
        std::size_t molecule_index;
        std::optional<std::size_t> conformer_index;
    };

    std::vector<Target> targets;

    std::size_t target_count = 0;
    for (const auto& prepared_molecule : prepared_collection.molecules()) {
        target_count += geometry_dependent ? prepared_molecule.molecule().conformer_count() : 1;
    }
    targets.reserve(target_count);

    for (std::size_t molecule_index = 0; molecule_index < prepared_collection.size();
         ++molecule_index) {
        const auto& prepared_molecule = prepared_collection[molecule_index];
        const auto& molecule = prepared_molecule.molecule();

        const auto molecule_target_count = geometry_dependent ? molecule.conformer_count() : 1;

        for (std::size_t target_index = 0; target_index < molecule_target_count; ++target_index) {
            targets.push_back(Target{
                .molecule_index = molecule_index,
                .conformer_index =
                    geometry_dependent ? std::optional<std::size_t>{target_index} : std::nullopt});
        }
    }

    std::vector<std::optional<charges::ChargeAssignment>> calculated(targets.size());
    ::chargefw::calculation::detail::parallel_for_indexed(
        targets.size(), max_threads, [&](const std::size_t target_index) {
            const auto& target = targets[target_index];
            const auto& prepared_molecule = prepared_collection[target.molecule_index];
            auto atomic_charges =
                selected.uses_parameters()
                    ? calculate_with_parameters(selected, prepared_molecule,
                                                selected.classifications[target.molecule_index],
                                                target.conformer_index)
                    : calculate_without_parameters(selected, prepared_molecule,
                                                   target.conformer_index);

            calculated[target_index] = charges::ChargeAssignment{
                .target = charges::ChargeTarget{.molecule_index = target.molecule_index,
                                                .conformer_index = target.conformer_index},
                .charges = std::move(atomic_charges)};
        });

    std::vector<charges::ChargeAssignment> assignments;
    assignments.reserve(calculated.size());
    for (auto& assignment : calculated) {
        assignments.push_back(std::move(*assignment));
    }

    return charges::ChargeSet{std::string{selected.method->id()}, std::move(assignments),
                              parameter_set_id_for(selected)};
}

} // namespace chargefw::methods
