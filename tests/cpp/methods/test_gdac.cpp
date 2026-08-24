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

#include <cmath>
#include <optional>
#include <snitch/snitch.hpp>
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

TEST_CASE("GDAC rejects missing features and responds to geometry", "[methods][gdac]") {
    const auto& registry = methods::method_registry();
    const auto* gdac = registry.find("gdac");

    REQUIRE(gdac != nullptr);
    const auto options = methods::make_default_options(gdac->option_schema());

    const auto parameter_set = make_test_gdac_parameters();

    const auto charged_pair = chargefw::test::make_formally_charged_pair();
    const features::PreparedMolecule prepared_charged_pair{charged_pair};

    const auto prerequisite_result = gdac->check_method_prerequisites(
        {.prepared_molecule = prepared_charged_pair, .method_options = options});

    CHECK(!prerequisite_result);
    REQUIRE(!prerequisite_result.issues().empty());
    CHECK(prerequisite_result.issues()[0].kind == methods::PrerequisiteIssueKind::missing_feature);

    const chargefw::core::Molecule rubidium_molecule{
        {chargefw::core::Atom{37}}, {}, {chargefw::core::Conformer{{chargefw::core::Position{}}}}};
    const features::PreparedMolecule prepared_rubidium{rubidium_molecule};
    const auto rubidium_prerequisite_result = gdac->check_method_prerequisites(
        {.prepared_molecule = prepared_rubidium, .method_options = options});

    CHECK(!rubidium_prerequisite_result);
    REQUIRE(rubidium_prerequisite_result.issues().size() == 1);
    CHECK(rubidium_prerequisite_result.issues()[0].kind ==
          methods::PrerequisiteIssueKind::unsupported_molecule);
    CHECK(rubidium_prerequisite_result.issues()[0].atom_index == 0);

    const auto two_conformer_water = chargefw::test::make_two_conformer_water();
    const auto two_conformer_charge_set =
        chargefw::test::calculate_method(two_conformer_water, "gdac", {parameter_set}, &options);

    const auto& first_assignment = two_conformer_charge_set.assignment(0);
    const auto& second_assignment = two_conformer_charge_set.assignment(1);
    const auto& first_charges = first_assignment.charges;
    const auto& second_charges = second_assignment.charges;

    CHECK(std::abs(first_charges[0] - second_charges[0]) > 1.0e-4);
}
