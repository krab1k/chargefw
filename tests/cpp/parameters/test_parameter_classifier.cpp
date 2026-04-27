#include "support/test_molecules.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/features/topology_features.h>
#include <chargefw/parameters/atom_parameters.h>
#include <chargefw/parameters/bond_parameters.h>
#include <chargefw/parameters/parameter_classifier.h>
#include <chargefw/parameters/parameter_key.h>
#include <chargefw/parameters/parameter_set.h>
#include <chargefw/parameters/parameter_set_metadata.h>

#include <cassert>
#include <string>
#include <utility>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace parameters = chargefw::parameters;

namespace {

auto test_metadata() -> parameters::ParameterSetMetadata {
    return {.id = "test-parameters", .method_id = "test-method", .name = "Test parameters"};
}

auto atom_key(const int atomic_number,
              const parameters::AtomParameterClassificationKind classification, std::string type)
    -> parameters::AtomParameterKey {
    return {
        .atomic_number = atomic_number, .classification = classification, .type = std::move(type)};
}

auto bond_type_key(const parameters::BondParameterClassificationKind classification,
                   std::string type) -> parameters::BondTypeKey {
    return {.classification = classification, .type = std::move(type)};
}

auto classify(const core::Molecule& molecule, const parameters::ParameterSet& parameter_set)
    -> parameters::ParameterClassification {
    const features::TopologyFeatures topology{molecule};

    return parameters::classify_parameters(molecule, topology, parameter_set);
}

auto make_c_cl_f_h_fragment() -> core::Molecule {
    std::vector atoms{core::Atom{6, 0, "C"}, core::Atom{17, 0, "Cl"}, core::Atom{9, 0, "F"},
                      core::Atom{1, 0, "H"}};

    std::vector bonds{core::Bond{0, 1, core::BondOrder::SINGLE},
                      core::Bond{0, 2, core::BondOrder::SINGLE},
                      core::Bond{0, 3, core::BondOrder::SINGLE}};

    return core::Molecule{std::move(atoms), std::move(bonds), {}, "c-cl-f-h"};
}

} // namespace

auto main() -> int {
    const auto water = chargefw::test::make_water();

    const parameters::ParameterSet water_parameters{
        test_metadata(),
        {},
        parameters::AtomParameters{
            {{.key = atom_key(1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key =
                  atom_key(8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}}}},
        parameters::BondParameters{
            {{.key = {.first_atom = atom_key(
                          1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
                      .second_atom = atom_key(
                          8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
                      .bond = bond_type_key(parameters::BondParameterClassificationKind::BOND_ORDER,
                                            "1")},
              .parameters = {{.name = "value", .value = 3.0}}}}}};

    const auto water_classification = classify(water, water_parameters);

    assert(water_classification.atom().size() == 3);
    assert(water_classification.atom()[0] == 1);
    assert(water_classification.atom()[1] == 0);
    assert(water_classification.atom()[2] == 0);

    assert(water_classification.bond().size() == 2);
    assert(water_classification.bond()[0] == 0);
    assert(water_classification.bond()[1] == 0);

    const auto mixed = make_c_cl_f_h_fragment();

    const parameters::ParameterSet mixed_parameters{
        test_metadata(),
        {},
        parameters::AtomParameters{
            {{.key =
                  atom_key(6, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "ClFH"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key =
                  atom_key(17, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "C"),
              .parameters = {{.name = "value", .value = 2.0}}},
             {.key = atom_key(9, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "C"),
              .parameters = {{.name = "value", .value = 3.0}}},
             {.key = atom_key(1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "C"),
              .parameters = {{.name = "value", .value = 4.0}}}}}};

    const auto mixed_classification = classify(mixed, mixed_parameters);

    assert(mixed_classification.atom().size() == 4);
    assert(mixed_classification.atom()[0] == 0);
    assert(mixed_classification.atom()[1] == 1);
    assert(mixed_classification.atom()[2] == 2);
    assert(mixed_classification.atom()[3] == 3);

    return 0;
}
