#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/core/molecule.h>
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

#include <cstddef>
#include <snitch/snitch.hpp>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

auto make_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-delre", .method_id = "delre", .name = "Test DelRe parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "delta", .value = 2.0}}},
                                    {.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "delta", .value = 1.0}}}}},
        parameters::BondParameters{{{.key = chargefw::test::plain_bond_key(8, 1),
                                     .parameters = {{.name = "eps", .value = 1.0},
                                                    {.name = "gammaA", .value = 0.2},
                                                    {.name = "gammaB", .value = 0.1}}}}}};
}

auto calculate_delre(const methods::Method& method, const core::Molecule& molecule,
                     const parameters::ParameterSet& parameter_set,
                     const parameters::ParameterClassification& classification)
    -> chargefw::charges::AtomicCharges {
    const features::PreparedMolecule prepared{molecule};
    const auto method_options = methods::make_default_options(method.option_schema());
    const parameters::ParameterView parameter_view{parameter_set, classification};
    const methods::CalculationInput input{
        prepared, method_options, core::total_formal_charge(molecule), nullptr, &parameter_view};

    return method.calculate(input);
}

} // namespace

TEST_CASE("DelRe produces conformer-independent water charges with explicit classification",
          "[methods][delre]") {
    const auto& registry = methods::method_registry();
    const auto* delre = registry.find("delre");

    CHECK(delre != nullptr);

    const auto parameter_set = make_parameter_set();

    const auto workflow_charge_set = chargefw::test::calculate_method(
        chargefw::test::make_water_graph(), "delre", {parameter_set});
    chargefw::test::assert_calculation_provenance(workflow_charge_set, "delre", "test-delre");
    chargefw::test::assert_conformer_independent(workflow_charge_set);

    {
        // Explicit classification pins parameter indices to source atom order; the charges
        // must follow that mapping exactly.
        const auto molecule = chargefw::test::make_water_graph();

        const auto classification = parameters::ParameterClassification{
            parameters::AtomParameterClassification{std::vector<std::size_t>{0, 1, 1}},
            parameters::BondParameterClassification{std::vector<std::size_t>{0, 0}}};

        const auto charges = calculate_delre(*delre, molecule, parameter_set, classification);

        CHECK(charges.size() == 3);
        CHECK(std::abs(charges[0] - (-1.25)) < 1.0e-12);
        CHECK(std::abs(charges[1] - (0.625)) < 1.0e-12);
        CHECK(std::abs(charges[2] - (0.625)) < 1.0e-12);
        CHECK(std::abs(charges.total() - (0.0)) < 1.0e-12);
    }

    chargefw::test::assert_water_charges_labeling_invariant("delre", {parameter_set});
    chargefw::test::assert_water_charges_geometry_independent("delre", {parameter_set});
}