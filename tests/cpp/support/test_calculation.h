#pragma once

#include "support/test_molecules.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

#include <snitch/snitch.hpp>

namespace chargefw::test {

[[nodiscard]] inline auto
calculate_method(core::Molecule molecule, std::string_view method_id,
                 std::vector<parameters::ParameterSet> parameter_sets = {},
                 const methods::MethodOptions* method_options = nullptr) -> charges::ChargeSet {
    const core::MoleculeCollection collection{std::vector{std::move(molecule)}, "test"};
    const features::PreparedMoleculeCollection prepared{collection};
    const auto* method = methods::method_registry().find(method_id);

    REQUIRE(method != nullptr);

    const std::vector<const methods::Method*> candidates{method};
    auto applicability = methods::find_applicable_methods(
        {.molecules = prepared, .methods = candidates, .parameter_sets = parameter_sets});

    REQUIRE(applicability.applicable.size() == 1);
    REQUIRE(applicability.rejected.empty());

    if (method_options != nullptr) {
        applicability.applicable.front().method_options = *method_options;
    }

    return calculation::calculate(
               {.molecules = prepared, .selected = applicability.applicable.front()})
        .charges;
}

[[nodiscard]] inline auto
calculate_single_method(core::Molecule molecule, std::string_view method_id,
                        std::vector<parameters::ParameterSet> parameter_sets = {},
                        const methods::MethodOptions* method_options = nullptr)
    -> charges::ChargeSet {
    const auto charge_set =
        calculate_method(std::move(molecule), method_id, std::move(parameter_sets), method_options);
    REQUIRE(charge_set.size() == 1);
    return charge_set;
}

inline auto assert_calculation_provenance(const charges::ChargeSet& charge_set,
                                          const std::string_view method_id,
                                          const std::optional<std::string_view> parameter_set_id)
    -> void {
    CHECK(charge_set.method_id() == method_id);
    CHECK(charge_set.parameter_set_id() == parameter_set_id);
}

inline auto assert_conformer_independent(const charges::ChargeSet& charge_set) -> void {
    CHECK(charge_set.size() == 1);
    const auto& assignment = charge_set.assignment(0);
    CHECK(assignment.target.molecule_index == 0);
    CHECK_FALSE(assignment.target.conformer_index.has_value());
}

inline auto assert_conformer_dependent(const charges::ChargeSet& charge_set,
                                       const std::size_t conformer_count) -> void {
    REQUIRE(charge_set.size() == conformer_count);

    for (std::size_t conformer_index = 0; conformer_index < conformer_count; ++conformer_index) {
        const auto& assignment = charge_set.assignment(conformer_index);
        CHECK(assignment.target.molecule_index == 0);
        CHECK(assignment.target.conformer_index == conformer_index);
    }
}

/// Asserts the chemically expected charge shape on neutral water: the oxygen carries the
/// negative charge, both hydrogens carry equal positive charge, and the molecule is neutral.
inline auto assert_neutral_water_charges(const charges::AtomicCharges& charges,
                                         const double tolerance) -> void {
    REQUIRE(charges.size() == 3);
    CHECK(charges[0] < 0.0);
    CHECK(charges[1] > 0.0);
    CHECK(charges[2] > 0.0);
    CHECK(std::abs(charges[1] - charges[2]) < tolerance);
    CHECK(std::abs(charges.total() - 0.0) < tolerance);
}

inline auto assert_same_charges(const charges::AtomicCharges& actual,
                                const charges::AtomicCharges& expected, const double tolerance)
    -> void {
    REQUIRE(actual.size() == expected.size());

    for (std::size_t atom_index = 0; atom_index < actual.size(); ++atom_index) {
        CHECK(std::abs(actual[atom_index] - expected[atom_index]) < tolerance);
    }
}

/// Method-agnostic invariant: charges must not depend on atom numbering or bond endpoint
/// order. Verifies a relabeled water (same physical molecule, permuted indices) and a
/// bond-flipped water against the plain-water reference.
inline auto
assert_water_charges_labeling_invariant(std::string_view method_id,
                                        std::vector<parameters::ParameterSet> parameter_sets = {},
                                        const methods::MethodOptions* method_options = nullptr,
                                        const double tolerance = 1.0e-9) -> void {
    static constexpr std::array new_to_old{std::size_t{2}, std::size_t{0}, std::size_t{1}};

    const auto reference =
        calculate_single_method(make_water(), method_id, parameter_sets, method_options);
    const auto& reference_charges = reference.assignment(0).charges;

    const auto relabeled = calculate_single_method(relabel_atoms(make_water(), new_to_old),
                                                   method_id, parameter_sets, method_options);
    const auto& relabeled_charges = relabeled.assignment(0).charges;

    REQUIRE(relabeled_charges.size() == reference_charges.size());
    for (std::size_t atom_index = 0; atom_index < new_to_old.size(); ++atom_index) {
        CHECK(std::abs(relabeled_charges[atom_index] -
                       reference_charges[new_to_old.at(atom_index)]) < tolerance);
    }

    const auto flipped = calculate_single_method(flip_bond_directions(make_water()), method_id,
                                                 std::move(parameter_sets), method_options);
    assert_same_charges(flipped.assignment(0).charges, reference_charges, tolerance);
}

/// Method-agnostic invariant for topology-only methods: geometry must not influence charges.
/// Verifies that a molecule with two distinct conformers yields a single conformer-independent
/// assignment identical to the conformer-free result.
inline auto
assert_water_charges_geometry_independent(std::string_view method_id,
                                          std::vector<parameters::ParameterSet> parameter_sets = {},
                                          const methods::MethodOptions* method_options = nullptr,
                                          const double tolerance = 1.0e-9) -> void {
    const auto reference =
        calculate_single_method(make_water_graph(), method_id, parameter_sets, method_options);

    const auto two_conformer = calculate_method(make_two_conformer_water(), method_id,
                                                std::move(parameter_sets), method_options);
    assert_conformer_independent(two_conformer);

    assert_same_charges(two_conformer.assignment(0).charges, reference.assignment(0).charges,
                        tolerance);
}

} // namespace chargefw::test
