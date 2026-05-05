#include "support/test_molecules.h"

#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/classification/parameter_classifier.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cassert>
#include <cmath>
#include <optional>
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

auto make_test_gdac_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "gdac-test", .method_id = "gdac", .name = "GDAC test parameters"},
        {},
        parameters::AtomParameters{
            {{.key = atom_key(1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "A", .value = 7.0}, {.name = "B", .value = 20.0}}},
             {.key =
                  atom_key(8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "A", .value = 12.0}, {.name = "B", .value = 10.0}}}}}};
}

} // namespace

auto main() -> int {
    const auto& registry = methods::method_registry();
    const auto* gdac = registry.find("gdac");

    assert(gdac != nullptr);
    assert(gdac->requires_parameters());
    assert(gdac->requirements().bond_graph);
    assert(gdac->requirements().coordinates);
    assert(gdac->requirements().element_properties);
    assert(gdac->requirements().atom_parameters.size() == 2);

    const auto options = methods::make_default_options(gdac->option_schema());

    const auto water = chargefw::test::make_water();
    const features::PreparedMolecule prepared_water{water};
    const features::ConformerFeatures geometry{water};

    const auto parameter_set = make_test_gdac_parameters();

    const auto classification =
        parameters::classify_parameters(water, prepared_water.topology(), parameter_set);

    const parameters::ParameterView parameter_view{parameter_set, classification};

    const methods::CalculationInput input{prepared_water, options, &geometry, &parameter_view};

    const auto charges = gdac->calculate(input);

    assert(charges.size() == water.atom_count());
    assert(charges[0] < 0.0);
    assert(charges[1] > 0.0);
    assert(charges[2] > 0.0);
    assert(std::abs(charges[1] - charges[2]) < 1.0e-4);
    assert(std::abs(charges.total()) < 1.0e-4);

    const auto charged_pair = chargefw::test::make_formally_charged_pair();
    const features::PreparedMolecule prepared_charged_pair{charged_pair};

    const auto prerequisite_result = gdac->check_method_prerequisites(
        {.prepared_molecule = prepared_charged_pair, .method_options = options});

    assert(!prerequisite_result);
    assert(!prerequisite_result.issues().empty());
    assert(prerequisite_result.issues()[0].kind == methods::PrerequisiteIssueKind::missing_feature);

    const auto two_conformer_water = chargefw::test::make_two_conformer_water();
    const features::PreparedMolecule prepared_two_conformer_water{two_conformer_water};

    const auto two_conformer_classification = parameters::classify_parameters(
        two_conformer_water, prepared_two_conformer_water.topology(), parameter_set);

    const parameters::ParameterView two_conformer_parameter_view{parameter_set,
                                                                 two_conformer_classification};

    const features::ConformerFeatures first_geometry{two_conformer_water, 0};
    const features::ConformerFeatures second_geometry{two_conformer_water, 1};

    const methods::CalculationInput first_input{prepared_two_conformer_water, options,
                                                &first_geometry, &two_conformer_parameter_view};

    const methods::CalculationInput second_input{prepared_two_conformer_water, options,
                                                 &second_geometry, &two_conformer_parameter_view};

    const auto first_charges = gdac->calculate(first_input);
    const auto second_charges = gdac->calculate(second_input);

    assert(first_charges.size() == two_conformer_water.atom_count());
    assert(second_charges.size() == two_conformer_water.atom_count());

    assert(std::abs(first_charges.total()) < 1.0e-4);
    assert(std::abs(second_charges.total()) < 1.0e-4);

    assert(std::abs(first_charges[0] - second_charges[0]) > 1.0e-4);

    return 0;
}