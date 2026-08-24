#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cstddef>
#include <snitch/snitch.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

auto make_water_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "test-water-parameters",
                                         .method_id = "test-method",
                                         .name = "Test water parameters"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}}}}};
}

} // namespace

TEST_CASE("calculation input exposes molecule, topology, geometry, and parameters",
          "[methods][calculation-input]") {
    const auto water = chargefw::test::make_water();
    const features::PreparedMolecule prepared_water{water};
    const methods::MethodOptions options;

    const methods::CalculationInput basic_input{prepared_water, options, -1.5};

    CHECK(&basic_input.prepared_molecule() == &prepared_water);
    CHECK(&basic_input.molecule() == &water);
    CHECK(&basic_input.topology() == &prepared_water.topology());
    CHECK(&basic_input.method_options() == &options);
    CHECK(basic_input.target_charge() == -1.5);

    CHECK_FALSE(basic_input.has_geometry());
    CHECK(basic_input.geometry_if_available() == nullptr);

    CHECK_THROWS_AS(basic_input.geometry(), std::logic_error);

    CHECK_FALSE(basic_input.has_parameters());
    CHECK(basic_input.parameters_if_available() == nullptr);

    CHECK_THROWS_AS(basic_input.parameters(), std::logic_error);

    const features::ConformerFeatures geometry{water};

    const methods::CalculationInput geometry_input{prepared_water, options, 0.0, &geometry};

    CHECK(geometry_input.has_geometry());
    CHECK(geometry_input.geometry_if_available() == &geometry);
    CHECK(&geometry_input.geometry() == &geometry);

    const auto parameter_set = make_water_parameters();

    const parameters::ParameterClassification classification{
        parameters::AtomParameterClassification{std::vector<std::size_t>{1, 0, 0}}};

    const parameters::ParameterView parameter_view{parameter_set, classification};

    const methods::CalculationInput parameter_input{prepared_water, options, 0.0, nullptr,
                                                    &parameter_view};

    CHECK(parameter_input.has_parameters());
    CHECK(parameter_input.parameters_if_available() == &parameter_view);
    CHECK(&parameter_input.parameters() == &parameter_view);

    const auto value = parameter_input.parameters().atom("value");

    CHECK(value[0] == 2.0);
    CHECK(value[1] == 1.0);
    CHECK(value[2] == 1.0);
}