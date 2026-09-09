#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
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

auto make_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-qeq", .method_id = "qeq", .name = "Test QEq parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 4.5280},
                                                    {.name = "hardness", .value = 13.8904}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 8.741},
                                                    {.name = "hardness", .value = 13.364}}}}}};
}

} // namespace

TEST_CASE("QEq defaults to DasGupta-Huzinaga and responds to geometry", "[methods][qeq]") {
    const auto charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "qeq", {make_parameter_set()});
    const auto& first_charges = charge_set.assignment(0).charges;
    const auto& second_charges = charge_set.assignment(1).charges;

    chargefw::test::assert_neutral_water_charges(first_charges, 1.0e-10);
    REQUIRE(second_charges.size() == 3);
    CHECK(second_charges[0] < 0.0);
    CHECK(second_charges[1] > 0.0);
    CHECK(second_charges[2] > 0.0);
    CHECK(std::abs(second_charges.total()) < 1.0e-10);
    CHECK(std::abs(first_charges[0] - second_charges[0]) > 1.0e-8);

    auto options = chargefw::methods::MethodOptions{};
    options.set("overlap_term", std::string{"DasGupta-Huzinaga"});
    const auto explicit_charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "qeq", {make_parameter_set()}, &options);

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
            chargefw::test::make_water(), "qeq", {make_parameter_set()}, &options);
        chargefw::test::assert_neutral_water_charges(charge_set.assignment(0).charges, 1.0e-10);
    }
}
