#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cmath>
#include <snitch/snitch.hpp>
#include <vector>

namespace parameters = chargefw::parameters;

namespace {

auto make_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-eem", .method_id = "eem", .name = "Test EEM parameters"},
        parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
        parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 5.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 9.0}}}}}};
}

} // namespace

TEST_CASE("EEM enforces a charged molecular target", "[methods][eem]") {
    const chargefw::core::Molecule cation{
        {chargefw::core::Atom{1, 1}},
        {},
        {chargefw::core::Conformer{{chargefw::core::Position{}}}}};
    const auto charge_set =
        chargefw::test::calculate_single_method(cation, "eem", {make_parameter_set()});

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
            chargefw::test::calculate_method(molecule, "eem", {make_parameter_set()});
        chargefw::test::assert_conformer_dependent(charge_set, 2);

        // EEM uses A = chi, B = hardness (no factor of two), and J = kappa/r.
        // Subtracting 5*q_H + J*q_O + lambda = -1 and
        // J*q_H + 9*q_O + lambda = -2 with q_H + q_O = Q gives
        // q_H = (1 + (9-J)*Q)/(14-2*J), q_O = Q-q_H.
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
