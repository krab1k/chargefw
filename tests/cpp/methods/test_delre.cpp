#include "support/test_assertions.h"
#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cassert>
#include <string_view>
#include <utility>
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
        parameters::AtomParameters{
            {{.key = {.atomic_number = 8,
                      .classification = parameters::AtomParameterClassificationKind::PLAIN,
                      .type = "*"},
              .parameters = {{.name = "delta", .value = 2.0}}},
             {.key = {.atomic_number = 1,
                      .classification = parameters::AtomParameterClassificationKind::PLAIN,
                      .type = "*"},
              .parameters = {{.name = "delta", .value = 1.0}}}}},
        parameters::BondParameters{
            {{.key = {.first_atom = {.atomic_number = 8,
                                     .classification =
                                         parameters::AtomParameterClassificationKind::PLAIN,
                                     .type = "*"},
                      .second_atom = {.atomic_number = 1,
                                      .classification =
                                          parameters::AtomParameterClassificationKind::PLAIN,
                                      .type = "*"},
                      .bond = {.classification = parameters::BondParameterClassificationKind::PLAIN,
                               .type = "*"}},
              .parameters = {{.name = "eps", .value = 1.0},
                             {.name = "gammaA", .value = 0.2},
                             {.name = "gammaB", .value = 0.1}}}}}};
}

auto make_water_oxygen_first() -> core::Molecule {
    std::vector atoms{core::Atom{8, 0, "O"}, core::Atom{1, 0, "H1"}, core::Atom{1, 0, "H2"}};

    std::vector bonds{core::Bond{0, 1, core::BondOrder::SINGLE},
                      core::Bond{0, 2, core::BondOrder::SINGLE}};

    return core::Molecule{std::move(atoms), std::move(bonds), {}, "water"};
}

auto make_water_hydrogen_first_bonds() -> core::Molecule {
    std::vector atoms{core::Atom{1, 0, "H1"}, core::Atom{8, 0, "O"}, core::Atom{1, 0, "H2"}};

    std::vector bonds{core::Bond{0, 1, core::BondOrder::SINGLE},
                      core::Bond{2, 1, core::BondOrder::SINGLE}};

    return core::Molecule{std::move(atoms), std::move(bonds), {}, "water-reversed"};
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

auto main() -> int {
    const auto& registry = methods::method_registry();
    const auto* delre = registry.find("delre");

    assert(delre != nullptr);

    const auto parameter_set = make_parameter_set();

    const auto workflow_charge_set = chargefw::test::calculate_method(
        chargefw::test::make_two_conformer_water(), "delre", {parameter_set});
    chargefw::test::assert_calculation_provenance(workflow_charge_set, "delre", "test-delre");
    chargefw::test::assert_conformer_independent(workflow_charge_set);

    {
        const auto molecule = make_water_oxygen_first();

        const auto classification = parameters::ParameterClassification{
            parameters::AtomParameterClassification{std::vector<std::size_t>{0, 1, 1}},
            parameters::BondParameterClassification{std::vector<std::size_t>{0, 0}}};

        const auto charges = calculate_delre(*delre, molecule, parameter_set, classification);

        assert(charges.size() == 3);
        chargefw::test::assert_close(charges[0], -1.25, 1.0e-12);
        chargefw::test::assert_close(charges[1], 0.625, 1.0e-12);
        chargefw::test::assert_close(charges[2], 0.625, 1.0e-12);
        chargefw::test::assert_close(charges.total(), 0.0, 1.0e-12);
    }

    {
        const auto molecule = make_water_hydrogen_first_bonds();

        const auto classification = parameters::ParameterClassification{
            parameters::AtomParameterClassification{std::vector<std::size_t>{1, 0, 1}},
            parameters::BondParameterClassification{std::vector<std::size_t>{0, 0}}};

        const auto charges = calculate_delre(*delre, molecule, parameter_set, classification);

        assert(charges.size() == 3);
        chargefw::test::assert_close(charges[0], 0.625, 1.0e-12);
        chargefw::test::assert_close(charges[1], -1.25, 1.0e-12);
        chargefw::test::assert_close(charges[2], 0.625, 1.0e-12);
        chargefw::test::assert_close(charges.total(), 0.0, 1.0e-12);
    }

    return 0;
}
