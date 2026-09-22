#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <array>
#include <cmath>
#include <snitch/snitch.hpp>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace parameters = chargefw::parameters;

namespace {

auto make_parameter_set(std::string_view method_id, std::string_view name)
    -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = std::string{"test-"} + std::string{method_id},
                                         .method_id = std::string{method_id},
                                         .name = std::string{name}},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 4.5280},
                                                    {.name = "hardness", .value = 13.8904},
                                                    {.name = "width", .value = 1.0}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 8.741},
                                                    {.name = "hardness", .value = 13.364},
                                                    {.name = "width", .value = 1.0}}}}},
        parameters::BondParameters{{{.key = chargefw::test::single_bond_key(1, 8),
                                     .parameters = {{.name = "kappa", .value = 1.0}}}}}};
}

} // namespace

TEST_CASE("SQE variants respond to changed conformer geometry", "[methods][sqe][sqeq0]") {
    constexpr auto variants = std::array{
        std::pair{"sqe", "Test SQE parameters"},
        std::pair{"sqeq0", "Test SQE+q0 parameters"},
    };

    for (const auto& [method_id, name] : variants) {
        CAPTURE(method_id);
        const auto charge_set =
            chargefw::test::calculate_method(chargefw::test::make_two_conformer_water(), method_id,
                                             {make_parameter_set(method_id, name)});
        const auto& charges = charge_set.assignment(0).charges;

        CHECK(std::abs(charges[0] - charge_set.assignment(1).charges[0]) > 1.0e-4);
    }
}
