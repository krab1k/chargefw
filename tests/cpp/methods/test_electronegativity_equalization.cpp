#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/methods/method_options.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <array>
#include <cmath>
#include <snitch/snitch.hpp>
#include <string>
#include <string_view>
#include <vector>

namespace parameters = chargefw::parameters;

namespace {

auto no_parameter_sets() -> std::vector<parameters::ParameterSet> {
    return {};
}

auto make_eem_parameters() -> std::vector<parameters::ParameterSet> {
    const auto parameter_set = parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-eem", .method_id = "eem", .name = "Test EEM parameters"},
        parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
        parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 5.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 9.0}}}}}};
    return {parameter_set};
}

auto make_qeq_parameters() -> std::vector<parameters::ParameterSet> {
    const auto parameter_set = parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-qeq", .method_id = "qeq", .name = "Test QEq parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 4.5280},
                                                    {.name = "hardness", .value = 13.8904}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 8.741},
                                                    {.name = "hardness", .value = 13.364}}}}}};
    return {parameter_set};
}

auto make_eqeqc_parameters() -> std::vector<parameters::ParameterSet> {
    const auto parameter_set = parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-eqeqc", .method_id = "eqeqc", .name = "Test EQeq+C parameters"},
        parameters::CommonParameters{{{.name = "alpha", .value = 1.0}}},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "Dz", .value = 0.1}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "Dz", .value = 0.2}}}}}};
    return {parameter_set};
}

auto make_sfkeem_parameters() -> std::vector<parameters::ParameterSet> {
    const auto parameter_set = parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-sfkeem", .method_id = "sfkeem", .name = "Test SFKEEM parameters"},
        parameters::CommonParameters{{{.name = "sigma", .value = 1.0}}},
        parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 10.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 10.0}}}}}};
    return {parameter_set};
}

auto make_abeem_parameters() -> std::vector<parameters::ParameterSet> {
    const auto parameter_set = parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-abeem", .method_id = "abeem", .name = "Test ABEEM parameters"},
        parameters::CommonParameters{{{.name = "k", .value = 1.0}}},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "a", .value = 1.0},
                                                    {.name = "b", .value = 10.0},
                                                    {.name = "c", .value = 0.5}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "a", .value = 2.0},
                                                    {.name = "b", .value = 10.0},
                                                    {.name = "c", .value = 0.5}}}}},
        parameters::BondParameters{{{.key = chargefw::test::plain_bond_key(8, 1),
                                     .parameters = {{.name = "A", .value = 1.0},
                                                    {.name = "B", .value = 10.0},
                                                    {.name = "C", .value = 0.5},
                                                    {.name = "D", .value = 0.5}}}}}};
    return {parameter_set};
}

auto make_smpqeq_parameters() -> std::vector<parameters::ParameterSet> {
    const auto parameter_set = parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-smpqeq", .method_id = "smpqeq", .name = "Test SMP/QEq parameters"},
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
    return {parameter_set};
}

struct GeometryMethodCase {
    std::string_view id;
    double minimum_change;
    std::vector<parameters::ParameterSet> (*make_parameters)();
};

} // namespace

TEST_CASE("electronegativity-equalization methods respond to changed conformer geometry",
          "[methods][eem][qeq][eqeq][eqeqc][sfkeem][abeem][smpqeq]") {
    constexpr auto methods = std::array{
        GeometryMethodCase{"eem", 1.0e-8, make_eem_parameters},
        GeometryMethodCase{"qeq", 1.0e-8, make_qeq_parameters},
        GeometryMethodCase{"eqeq", 1.0e-8, no_parameter_sets},
        GeometryMethodCase{"eqeqc", 1.0e-8, make_eqeqc_parameters},
        GeometryMethodCase{"sfkeem", 1.0e-8, make_sfkeem_parameters},
        GeometryMethodCase{"abeem", 1.0e-4, make_abeem_parameters},
        GeometryMethodCase{"smpqeq", 1.0e-8, make_smpqeq_parameters},
    };

    for (const auto& method : methods) {
        CAPTURE(method.id);
        const auto charge_set = chargefw::test::calculate_method(
            chargefw::test::make_two_conformer_water(), method.id, method.make_parameters());

        CHECK(std::abs(charge_set.assignment(0).charges[0] - charge_set.assignment(1).charges[0]) >
              method.minimum_change);
    }
}

TEST_CASE("EEM enforces a charged molecular target", "[methods][eem]") {
    const chargefw::core::Molecule cation{
        {chargefw::core::Atom{1, 1}},
        {},
        {chargefw::core::Conformer{{chargefw::core::Position{}}}}};
    const auto charge_set =
        chargefw::test::calculate_single_method(cation, "eem", make_eem_parameters());

    CHECK(std::abs(charge_set.assignment(0).charges[0] - 1.0) < 1.0e-12);
    CHECK(std::abs(charge_set.assignment(0).charges.total() - 1.0) < 1.0e-12);
}

TEST_CASE("EEM two-atom charges obey equalization at different distances and targets",
          "[methods][eem]") {
    for (const int target : {0, 1}) {
        CAPTURE(target);
        const chargefw::core::Molecule molecule{
            {chargefw::core::Atom{1, target}, chargefw::core::Atom{8}},
            {},
            {chargefw::core::Conformer{
                 {chargefw::core::Position{}, chargefw::core::Position{.x = 1.0}}},
             chargefw::core::Conformer{
                 {chargefw::core::Position{}, chargefw::core::Position{.x = 2.0}}}}};
        const auto charge_set =
            chargefw::test::calculate_method(molecule, "eem", make_eem_parameters());
        chargefw::test::assert_conformer_dependent(charge_set, 2);

        // EEM uses A = chi, B = hardness (no factor of two), and J = kappa/r.
        // The equalization equations give q_H = (1 + (9-J)*Q)/(14-2*J) and q_O = Q-q_H.
        // At r=1: Q=0 -> (1/12,-1/12), Q=1 -> (3/4,1/4).
        // At r=2: Q=0 -> (1/13,-1/13), Q=1 -> (19/26,7/26).
        const std::vector expected_hydrogen{target == 0 ? 1.0 / 12.0 : 3.0 / 4.0,
                                            target == 0 ? 1.0 / 13.0 : 19.0 / 26.0};
        for (std::size_t index = 0; index < expected_hydrogen.size(); ++index) {
            CAPTURE(index);
            const auto& charges = charge_set.assignment(index).charges;
            REQUIRE(charges.size() == 2);
            CHECK(std::abs(charges[0] - expected_hydrogen[index]) < 1.0e-12);
            CHECK(std::abs(charges[1] - (target - expected_hydrogen[index])) < 1.0e-12);
            CHECK(std::abs(charges.total() - target) < 1.0e-12);
        }
    }
}

TEST_CASE("QEq defaults to DasGupta-Huzinaga", "[methods][qeq]") {
    const auto charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "qeq", make_qeq_parameters());
    const auto& first_charges = charge_set.assignment(0).charges;
    const auto& second_charges = charge_set.assignment(1).charges;

    chargefw::test::assert_neutral_water_charges(first_charges, 1.0e-10);
    REQUIRE(second_charges.size() == 3);
    CHECK(second_charges[0] < 0.0);
    CHECK(second_charges[1] > 0.0);
    CHECK(second_charges[2] > 0.0);
    CHECK(std::abs(second_charges.total()) < 1.0e-10);

    auto options = chargefw::methods::MethodOptions{};
    options.set("overlap_term", std::string{"DasGupta-Huzinaga"});
    const auto explicit_charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "qeq", make_qeq_parameters(), &options);

    chargefw::test::assert_same_charges(explicit_charge_set.assignment(0).charges, first_charges,
                                        1.0e-12);
    chargefw::test::assert_same_charges(explicit_charge_set.assignment(1).charges, second_charges,
                                        1.0e-12);
}

TEST_CASE("QEq empirical Coulomb terms calculate neutral water", "[methods][qeq]") {
    static constexpr std::array terms{
        std::string_view{"Nishimoto-Mataga"},
        std::string_view{"Nishimoto-Mataga-Weiss"},
        std::string_view{"Ohno"},
        std::string_view{"Ohno-Klopman"},
        std::string_view{"DasGupta-Huzinaga"},
        std::string_view{"Louwen-Vogt"},
    };

    for (const auto term : terms) {
        CAPTURE(term);
        auto options = chargefw::methods::MethodOptions{};
        options.set("overlap_term", std::string{term});

        const auto charge_set = chargefw::test::calculate_single_method(
            chargefw::test::make_water(), "qeq", make_qeq_parameters(), &options);
        chargefw::test::assert_neutral_water_charges(charge_set.assignment(0).charges, 1.0e-10);
    }
}
