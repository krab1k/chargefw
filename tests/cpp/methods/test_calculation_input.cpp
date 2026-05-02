#include "support/test_molecules.h"

#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/parameters/atom_parameters.h>
#include <chargefw/parameters/parameter_classification.h>
#include <chargefw/parameters/parameter_key.h>
#include <chargefw/parameters/parameter_set.h>
#include <chargefw/parameters/parameter_set_metadata.h>
#include <chargefw/parameters/parameter_view.h>

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

auto atom_key(const int atomic_number,
              const parameters::AtomParameterClassificationKind classification, std::string type)
    -> parameters::AtomParameterKey {
    return {
        .atomic_number = atomic_number, .classification = classification, .type = std::move(type)};
}

auto make_water_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "test-water-parameters",
                                         .method_id = "test-method",
                                         .name = "Test water parameters"},
        {},
        parameters::AtomParameters{
            {{.key = atom_key(1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key =
                  atom_key(8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}}}}};
}

} // namespace

auto main() -> int {
    const auto water = chargefw::test::make_water();
    const features::PreparedMolecule prepared_water{water};
    const methods::MethodOptions options;

    const methods::CalculationInput basic_input{prepared_water, options};

    assert(&basic_input.prepared_molecule() == &prepared_water);
    assert(&basic_input.molecule() == &water);
    assert(&basic_input.topology() == &prepared_water.topology());
    assert(&basic_input.method_options() == &options);

    assert(!basic_input.has_geometry());
    assert(basic_input.geometry_if_available() == nullptr);

    bool rejected_missing_geometry = false;

    try {
        [[maybe_unused]] const auto& geometry = basic_input.geometry();
    } catch (const std::logic_error&) {
        rejected_missing_geometry = true;
    }

    assert(rejected_missing_geometry);

    assert(!basic_input.has_parameters());
    assert(basic_input.parameters_if_available() == nullptr);

    bool rejected_missing_parameters = false;

    try {
        [[maybe_unused]] const auto& parameters = basic_input.parameters();
    } catch (const std::logic_error&) {
        rejected_missing_parameters = true;
    }

    assert(rejected_missing_parameters);

    const features::ConformerFeatures geometry{water};

    const methods::CalculationInput geometry_input{prepared_water, options, &geometry};

    assert(geometry_input.has_geometry());
    assert(geometry_input.geometry_if_available() == &geometry);
    assert(&geometry_input.geometry() == &geometry);

    const auto parameter_set = make_water_parameters();

    const parameters::ParameterClassification classification{
        parameters::AtomParameterClassification{std::vector<std::size_t>{1, 0, 0}}};

    const parameters::ParameterView parameter_view{parameter_set, classification};

    const methods::CalculationInput parameter_input{prepared_water, options, nullptr,
                                                    &parameter_view};

    assert(parameter_input.has_parameters());
    assert(parameter_input.parameters_if_available() == &parameter_view);
    assert(&parameter_input.parameters() == &parameter_view);

    const auto value = parameter_input.parameters().atom("value");

    assert(value[0] == 2.0);
    assert(value[1] == 1.0);
    assert(value[2] == 1.0);

    return 0;
}