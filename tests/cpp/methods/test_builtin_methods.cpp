#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <array>
#include <cmath>
#include <snitch/snitch.hpp>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

struct MethodManifest {
    std::string_view id;
    std::string_view name;
    std::string_view full_name;
    bool has_publication;
    int priority;
    bool coordinates;
    std::size_t common_parameter_count;
    std::size_t atom_parameter_count;
    std::size_t bond_parameter_count;
    std::size_t option_count;
    methods::ComplexityTerm time;
    methods::ComplexityTerm memory;
    bool supports_cutoff;
    bool supports_cover;
    methods::FragmentTargetChargePolicy fragment_target_charge_policy;
    bool uses_bundled_parameters = true;
};

using Complexity = methods::ComplexityTerm;
using FragmentCharge = methods::FragmentTargetChargePolicy;

#ifndef CHARGEFW_TEST_PARAMETER_DIR
#error "CHARGEFW_TEST_PARAMETER_DIR must be defined"
#endif

constexpr std::array method_manifest{
    MethodManifest{"abeem", "ABEEM", "Atom-Bond Electronegativity Equalization Method", true, 190,
                   true, 1, 3, 4, 0, Complexity::atoms_plus_bonds_cubed,
                   Complexity::atoms_plus_bonds_squared, true, true,
                   FragmentCharge::proportional_to_atom_count},
    MethodManifest{"charge2", "Charge2", "Charge2", true, 30, false, 6, 3, 0, 1,
                   Complexity::atoms_plus_bonds, Complexity::constant, false, false,
                   FragmentCharge::unsupported},
    MethodManifest{"delre", "DelRe", "Method of Del Re", true, 130, false, 0, 1, 3, 0,
                   Complexity::atoms_cubed, Complexity::atoms_squared, false, false,
                   FragmentCharge::unsupported},
    MethodManifest{"denr", "DENR", "Dynamical Electronegativity Relaxation", true, 50, false, 0, 2,
                   0, 2, Complexity::atoms_cubed, Complexity::atoms_squared, false, false,
                   FragmentCharge::unsupported},
    MethodManifest{"dummy", "Dummy method", "Dummy zero charges", false, 0, false, 0, 0, 0, 0,
                   Complexity::atoms, Complexity::constant, false, false,
                   FragmentCharge::unsupported},
    MethodManifest{"eem", "EEM", "Electronegativity Equalization Method", true, 200, true, 1, 2, 0,
                   0, Complexity::atoms_cubed, Complexity::atoms_squared, true, true,
                   FragmentCharge::proportional_to_atom_count},
    MethodManifest{"eqeq", "EQeq", "Extended Charge Equilibration Method", true, 150, true, 0, 0, 0,
                   0, Complexity::atoms_cubed, Complexity::atoms_squared, true, true,
                   FragmentCharge::proportional_to_atom_count},
    MethodManifest{"eqeqc", "EQeq+C", "Bond-Order-Corrected Extended Charge Equilibration Method",
                   true, 140, true, 1, 1, 0, 0, Complexity::atoms_cubed, Complexity::atoms_squared,
                   true, true, FragmentCharge::proportional_to_atom_count},
    MethodManifest{"formal", "Formal", "Formal atomic charges", false, 10, false, 0, 0, 0, 0,
                   Complexity::atoms, Complexity::constant, false, false,
                   FragmentCharge::unsupported},
    MethodManifest{"gdac", "GDAC", "Geometry-Dependent Net Atomic Charges", true, 100, true, 0, 2,
                   0, 1, Complexity::atoms_plus_bonds, Complexity::atoms, false, false,
                   FragmentCharge::unsupported},
    MethodManifest{"kcm", "KCM", "Kirchhoff Charge Model", true, 60, false, 0, 2, 0, 0,
                   Complexity::atoms_cubed, Complexity::atoms_squared, false, false,
                   FragmentCharge::unsupported},
    MethodManifest{"mgc", "MGC", "Molecular Graph Charge", true, 70, false, 0, 0, 0, 0,
                   Complexity::atoms_cubed, Complexity::atoms_squared, false, false,
                   FragmentCharge::unsupported},
    MethodManifest{"mpeoe", "MPEOE", "Modified Partial Equalization of Atomic Electronegativity",
                   true, 110, false, 1, 2, 1, 1, Complexity::atoms_plus_bonds, Complexity::atoms,
                   false, false, FragmentCharge::unsupported},
    MethodManifest{"peoe", "PEOE", "Partial Equalization of Atomic Electronegativity", true, 120,
                   false, 1, 3, 0, 1, Complexity::atoms_plus_bonds, Complexity::atoms, false, false,
                   FragmentCharge::unsupported},
    MethodManifest{"qeq", "QEq", "Charge Equilibration", true, 170, true, 0, 2, 0, 1,
                   Complexity::atoms_cubed, Complexity::atoms_squared, true, true,
                   FragmentCharge::proportional_to_atom_count},
    MethodManifest{"sfkeem", "SFKEEM",
                   "Selfconsistent Functional Kernel Equalized Electronegativity Method", true, 180,
                   true, 1, 2, 0, 0, Complexity::atoms_cubed, Complexity::atoms_squared, true, true,
                   FragmentCharge::proportional_to_atom_count},
    MethodManifest{"smpqeq", "SMP/QEq", "Self-Consistent Charge Equilibration Method", true, 160,
                   true, 0, 4, 0, 0, Complexity::atoms_cubed, Complexity::atoms_squared, false,
                   false, FragmentCharge::unsupported, false},
    MethodManifest{"sqe", "SQE", "Split-charge Equilibration", true, 90, true, 0, 3, 1, 0,
                   Complexity::bonds_cubed, Complexity::bonds_squared, true, true,
                   FragmentCharge::zero},
    MethodManifest{"sqeq0", "SQE+q0", "Split-charge Equilibration with Initial Formal Charges",
                   true, 80, true, 0, 3, 1, 0, Complexity::bonds_cubed, Complexity::bonds_squared,
                   true, true, FragmentCharge::proportional_to_atom_count},
    MethodManifest{"sqeqp", "SQE+qp",
                   "Split-charge Equilibration with Parameterized Initial Charges", true, 210, true,
                   0, 4, 1, 0, Complexity::bonds_cubed, Complexity::bonds_squared, true, true,
                   FragmentCharge::proportional_to_atom_count},
    MethodManifest{"tsef", "TSEF", "Topologically Symmetrical Energy Function", true, 55, false, 0,
                   2, 0, 0, Complexity::atoms_cubed, Complexity::atoms_squared, false, false,
                   FragmentCharge::unsupported},
    MethodManifest{"veem", "VEEM", "Valence Electrons Equalization Method", true, 20, false, 0, 0,
                   0, 0, Complexity::atoms, Complexity::constant, false, false,
                   FragmentCharge::unsupported},
};

auto calculate(const methods::Method& method, const chargefw::core::Molecule& molecule)
    -> chargefw::charges::AtomicCharges {
    const features::PreparedMolecule prepared_molecule{molecule};
    const auto method_options = methods::make_default_options(method.option_schema());

    const methods::CalculationInput input{prepared_molecule, method_options,
                                          chargefw::core::total_formal_charge(molecule)};

    return method.calculate(input);
}

auto make_smpqeq_water_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-smpqeq", .method_id = "smpqeq", .name = "Test SMP/QEq water parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "first", .value = 1.0},
                                                    {.name = "second", .value = 10.0},
                                                    {.name = "third", .value = 1.0},
                                                    {.name = "fourth", .value = 0.1}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "first", .value = 2.0},
                                                    {.name = "second", .value = 10.0},
                                                    {.name = "third", .value = 1.0},
                                                    {.name = "fourth", .value = 0.1}}}}}};
}

} // namespace

TEST_CASE("built-in method registry matches the conformance manifest",
          "[methods][builtin-methods]") {
    const auto& registry = methods::method_registry();
    const auto method_names = registry.names();

    REQUIRE(method_names.size() == method_manifest.size());

    for (const auto& expected : method_manifest) {
        const auto* method = registry.find(expected.id);
        REQUIRE(method != nullptr);

        CHECK(method->id() == expected.id);
        CHECK(method->metadata().name == expected.name);
        CHECK(method->metadata().full_name == expected.full_name);
        CHECK(method->metadata().publication.has_value() == expected.has_publication);
        CHECK(method->metadata().priority == expected.priority);

        const auto requirements = method->requirements();
        CHECK(requirements.coordinates == expected.coordinates);
        CHECK(requirements.common_parameters.size() == expected.common_parameter_count);
        CHECK(requirements.atom_parameters.size() == expected.atom_parameter_count);
        CHECK(requirements.bond_parameters.size() == expected.bond_parameter_count);
        CHECK(method->requires_parameters() ==
              (expected.common_parameter_count + expected.atom_parameter_count +
                   expected.bond_parameter_count >
               0));
        CHECK(method->option_schema().size() == expected.option_count);
        CHECK(requirements.resources.time == expected.time);
        CHECK(requirements.resources.memory == expected.memory);
        CHECK(requirements.resources.supports_cutoff == expected.supports_cutoff);
        CHECK(requirements.resources.supports_cover == expected.supports_cover);
        CHECK(requirements.resources.fragment_target_charge_policy ==
              expected.fragment_target_charge_policy);
    }
}

TEST_CASE("DENR exposes its step and iteration defaults", "[methods][builtin-methods]") {
    const auto* denr = methods::method_registry().find("denr");
    REQUIRE(denr != nullptr);

    const auto schema = denr->option_schema();
    REQUIRE(schema.size() == 2);
    CHECK(schema[0].id == std::string_view{"step"});
    CHECK(schema[0].type == methods::MethodOptionType::floating_point);
    REQUIRE(std::get_if<double>(&schema[0].default_value) != nullptr);
    CHECK(*std::get_if<double>(&schema[0].default_value) == 0.1);
    CHECK(schema[1].id == std::string_view{"iterations"});
    CHECK(schema[1].type == methods::MethodOptionType::integer);
    REQUIRE(std::get_if<int>(&schema[1].default_value) != nullptr);
    CHECK(*std::get_if<int>(&schema[1].default_value) == 3);
}

TEST_CASE("every built-in completes its declared full workflow", "[methods][builtin-methods]") {
    const auto molecule = chargefw::test::make_two_conformer_water();
    const auto collection = chargefw::core::MoleculeCollection{std::vector{molecule}, "water"};
    const auto prepared = features::PreparedMoleculeCollection{collection};
    const auto parameter_sets =
        parameters::load_parameter_sets_json_directory(CHARGEFW_TEST_PARAMETER_DIR);
    const std::array smpqeq_parameter_sets{make_smpqeq_water_parameters()};

    for (const auto& expected : method_manifest) {
        CAPTURE(expected.id);
        const auto* method = methods::method_registry().find(expected.id);
        REQUIRE(method != nullptr);

        const std::array candidates{method};
        const auto candidate_parameter_sets =
            expected.uses_bundled_parameters
                ? std::span<const parameters::ParameterSet>{parameter_sets}
                : std::span<const parameters::ParameterSet>{smpqeq_parameter_sets};
        const auto applicability =
            methods::find_applicable_methods({.molecules = prepared,
                                              .methods = candidates,
                                              .parameter_sets = candidate_parameter_sets});
        REQUIRE_FALSE(applicability.applicable.empty());
        const auto& selected = applicability.applicable.front();

        const auto result =
            chargefw::calculation::calculate({.molecules = prepared, .selected = selected});
        CHECK(result.charges.method_id() == expected.id);
        CHECK(result.charges.parameter_set_id().has_value() == method->requires_parameters());
        CHECK(result.charges.size() == (expected.coordinates ? molecule.conformer_count() : 1));

        for (const auto& assignment : result.charges.assignments()) {
            CHECK(assignment.target.molecule_index == 0);
            CHECK(assignment.charges.size() == molecule.atom_count());
            for (const auto charge : assignment.charges.values()) {
                CHECK(std::isfinite(charge));
            }
            CHECK(std::abs(assignment.charges.total()) < 1.0e-8);
        }

        if (!expected.coordinates) {
            const auto graph_water = chargefw::test::make_water_graph();
            const auto graph_collection =
                chargefw::core::MoleculeCollection{std::vector{graph_water}, "water"};
            const auto graph_prepared = features::PreparedMoleculeCollection{graph_collection};
            const auto graph_applicability =
                methods::find_applicable_methods({.molecules = graph_prepared,
                                                  .methods = candidates,
                                                  .parameter_sets = candidate_parameter_sets});
            REQUIRE_FALSE(graph_applicability.applicable.empty());
            const auto& graph_selected = graph_applicability.applicable.front();

            const auto graph_result = chargefw::calculation::calculate(
                {.molecules = graph_prepared, .selected = graph_selected});
            REQUIRE(graph_result.charges.size() == 1);
            chargefw::test::assert_same_charges(graph_result.charges.assignment(0).charges,
                                                result.charges.assignment(0).charges, 1.0e-10);
        }

        if (expected.id != "formal" && expected.id != "dummy") {
            const auto water = chargefw::test::make_water();
            const auto water_collection =
                chargefw::core::MoleculeCollection{std::vector{water}, "water"};
            const auto water_prepared = features::PreparedMoleculeCollection{water_collection};
            const auto water_applicability =
                methods::find_applicable_methods({.molecules = water_prepared,
                                                  .methods = candidates,
                                                  .parameter_sets = candidate_parameter_sets});
            REQUIRE_FALSE(water_applicability.applicable.empty());
            const auto& water_selected = water_applicability.applicable.front();

            const auto water_result = chargefw::calculation::calculate(
                {.molecules = water_prepared, .selected = water_selected});
            REQUIRE(water_result.charges.size() == 1);
            chargefw::test::assert_neutral_water_charges(water_result.charges.assignment(0).charges,
                                                         1.0e-10);
        }

        const std::vector<parameters::ParameterSet> selected_parameter_sets =
            selected.parameter_set == nullptr
                ? std::vector<parameters::ParameterSet>{}
                : std::vector<parameters::ParameterSet>{*selected.parameter_set};
        chargefw::test::assert_water_charges_labeling_invariant(expected.id,
                                                                selected_parameter_sets);
    }
}

TEST_CASE("dummy assigns exact zero charges", "[methods][builtin-methods]") {
    const auto* dummy = methods::method_registry().find("dummy");
    REQUIRE(dummy != nullptr);

    const auto water = chargefw::test::make_water_graph();
    const auto dummy_charges = calculate(*dummy, water);
    CHECK(dummy_charges.size() == water.atom_count());

    for (const auto charge : dummy_charges.values()) {
        CHECK(charge == 0.0);
    }
}

TEST_CASE("formal copies atomic formal charges", "[methods][builtin-methods]") {
    const auto* formal = methods::method_registry().find("formal");
    REQUIRE(formal != nullptr);

    const auto charged_pair = chargefw::test::make_formally_charged_pair();
    const auto formal_charges = calculate(*formal, charged_pair);
    REQUIRE(formal_charges.size() == charged_pair.atom_count());
    CHECK(formal_charges[0] == 1.0);
    CHECK(formal_charges[1] == -1.0);
    CHECK(formal_charges.total() == 0.0);
}
