#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cmath>
#include <snitch/snitch.hpp>
#include <vector>

namespace parameters = chargefw::parameters;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

auto make_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-sqeqp", .method_id = "sqeqp", .name = "Test SQE+qp parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 4.5280},
                                                    {.name = "hardness", .value = 13.8904},
                                                    {.name = "width", .value = 1.0},
                                                    {.name = "q0", .value = 0.25}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 8.741},
                                                    {.name = "hardness", .value = 13.364},
                                                    {.name = "width", .value = 1.0},
                                                    {.name = "q0", .value = -0.5}}}}},
        parameters::BondParameters{{{.key = chargefw::test::single_bond_key(1, 8),
                                     .parameters = {{.name = "kappa", .value = 1.0}}}}}};
}

} // namespace

TEST_CASE("SQE+qp responds to changed conformer geometry", "[methods][sqeqp]") {
    const auto charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "sqeqp", {make_parameter_set()});
    const auto& charges = charge_set.assignment(0).charges;

    CHECK(std::abs(charges[0] - charge_set.assignment(1).charges[0]) > 1.0e-4);
}

TEST_CASE("SQE+qp normalizes initial charges to the calculation target", "[methods][sqeqp]") {
    const auto molecule = chargefw::test::make_water();
    const features::PreparedMolecule prepared{molecule};
    const features::ConformerFeatures geometry{molecule, 0};
    const auto parameter_set = make_parameter_set();
    const auto classification = parameters::ParameterClassification{
        parameters::AtomParameterClassification{std::vector<std::size_t>{1, 0, 0}},
        parameters::BondParameterClassification{std::vector<std::size_t>{0, 0}}};
    const parameters::ParameterView parameter_view{parameter_set, classification};
    const auto* method = methods::method_registry().find("sqeqp");
    REQUIRE(method != nullptr);
    const auto options = methods::make_default_options(method->option_schema());
    constexpr auto target_charge = 1.25;
    const methods::CalculationInput input{prepared, options, target_charge, &geometry,
                                          &parameter_view};

    const auto charges = method->calculate(input);

    CHECK(std::abs(charges.total() - target_charge) < 1.0e-10);
}
