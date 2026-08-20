#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
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

auto make_test_gdac_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "gdac-test", .method_id = "gdac", .name = "GDAC test parameters"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "A", .value = 7.0}, {.name = "B", .value = 20.0}}},
             {.key = chargefw::test::atom_key(
                  8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "A", .value = 12.0}, {.name = "B", .value = 10.0}}}}}};
}

} // namespace

auto main() -> int {
    const auto& registry = methods::method_registry();
    const auto* gdac = registry.find("gdac");

    assert(gdac != nullptr);
    assert(gdac->requires_parameters());
    assert(gdac->requirements().coordinates);
    assert(gdac->requirements().atom_parameters.size() == 2);

    const auto options = methods::make_default_options(gdac->option_schema());

    const auto parameter_set = make_test_gdac_parameters();
    const auto charge_set = chargefw::test::calculate_single_method(
        chargefw::test::make_water(), "gdac", {parameter_set}, &options);
    const auto& charges = charge_set.assignment(0).charges;

    chargefw::test::assert_calculation_provenance(charge_set, "gdac", "gdac-test");
    assert(charges.size() == 3);
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
    const auto two_conformer_charge_set =
        chargefw::test::calculate_method(two_conformer_water, "gdac", {parameter_set}, &options);

    chargefw::test::assert_calculation_provenance(two_conformer_charge_set, "gdac", "gdac-test");
    chargefw::test::assert_conformer_dependent(two_conformer_charge_set, 2);

    const auto& first_assignment = two_conformer_charge_set.assignment(0);
    const auto& second_assignment = two_conformer_charge_set.assignment(1);
    const auto& first_charges = first_assignment.charges;
    const auto& second_charges = second_assignment.charges;

    assert(first_charges.size() == two_conformer_water.atom_count());
    assert(second_charges.size() == two_conformer_water.atom_count());

    assert(std::abs(first_charges.total()) < 1.0e-4);
    assert(std::abs(second_charges.total()) < 1.0e-4);

    assert(std::abs(first_charges[0] - second_charges[0]) > 1.0e-4);

    return 0;
}
