#pragma once

#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <cassert>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::test {

[[nodiscard]] inline auto
calculate_method(core::Molecule molecule, std::string_view method_id,
                 std::vector<parameters::ParameterSet> parameter_sets = {},
                 const methods::MethodOptions* method_options = nullptr) -> charges::ChargeSet {
    const core::MoleculeCollection collection{std::vector{std::move(molecule)}, "test"};
    const features::PreparedMoleculeCollection prepared{collection};
    const auto* method = methods::method_registry().find(method_id);

    assert(method != nullptr);

    const std::vector<const methods::Method*> candidates{method};
    auto applicability = methods::find_applicable_methods(
        {.molecules = prepared, .methods = candidates, .parameter_sets = parameter_sets});

    assert(applicability.applicable.size() == 1);
    assert(applicability.rejected.empty());

    if (method_options != nullptr) {
        applicability.applicable.front().method_options = *method_options;
    }

    return methods::calculate_charges(applicability.applicable.front(), prepared);
}

[[nodiscard]] inline auto
calculate_single_method(core::Molecule molecule, std::string_view method_id,
                        std::vector<parameters::ParameterSet> parameter_sets = {},
                        const methods::MethodOptions* method_options = nullptr)
    -> charges::ChargeSet {
    const auto charge_set =
        calculate_method(std::move(molecule), method_id, std::move(parameter_sets), method_options);
    assert(charge_set.size() == 1);
    return charge_set;
}

inline auto assert_calculation_provenance(const charges::ChargeSet& charge_set,
                                          const std::string_view method_id,
                                          const std::optional<std::string_view> parameter_set_id)
    -> void {
    assert(charge_set.method_id() == method_id);
    assert(charge_set.parameter_set_id() == parameter_set_id);
}

inline auto assert_conformer_independent(const charges::ChargeSet& charge_set) -> void {
    assert(charge_set.size() == 1);
    const auto& assignment = charge_set.assignment(0);
    assert(assignment.target.molecule_index == 0);
    assert(!assignment.target.conformer_index.has_value());
}

inline auto assert_conformer_dependent(const charges::ChargeSet& charge_set,
                                       const std::size_t conformer_count) -> void {
    assert(charge_set.size() == conformer_count);

    for (std::size_t conformer_index = 0; conformer_index < conformer_count; ++conformer_index) {
        const auto& assignment = charge_set.assignment(conformer_index);
        assert(assignment.target.molecule_index == 0);
        assert(assignment.target.conformer_index == conformer_index);
    }
}

} // namespace chargefw::test
