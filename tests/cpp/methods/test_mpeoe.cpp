#include "support/test_calculation.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <cassert>
#include <cmath>
#include <filesystem>
#include <utility>
#include <vector>

namespace core = chargefw::core;
namespace parameters = chargefw::parameters;

namespace {

#ifndef CHARGEFW_TEST_PARAMETER_DIR
#error "CHARGEFW_TEST_PARAMETER_DIR must be defined"
#endif

auto make_methane() -> core::Molecule {
    std::vector atoms{core::Atom{6, 0, "C"}, core::Atom{1, 0, "H1"}, core::Atom{1, 0, "H2"},
                      core::Atom{1, 0, "H3"}, core::Atom{1, 0, "H4"}};

    std::vector bonds{
        core::Bond{0, 1, core::BondOrder::SINGLE}, core::Bond{0, 2, core::BondOrder::SINGLE},
        core::Bond{0, 3, core::BondOrder::SINGLE}, core::Bond{0, 4, core::BondOrder::SINGLE}};

    return core::Molecule{std::move(atoms), std::move(bonds), {}, "methane"};
}

auto load_mpeoe_parameters() -> std::vector<parameters::ParameterSet> {
    const auto path = std::filesystem::path{CHARGEFW_TEST_PARAMETER_DIR} / "MPEOE_original.json";

    std::vector<parameters::ParameterSet> parameter_sets;
    parameter_sets.push_back(parameters::load_parameter_set_json_file(path));
    return parameter_sets;
}

} // namespace

auto main() -> int {
    const auto charge_set =
        chargefw::test::calculate_single_method(make_methane(), "mpeoe", load_mpeoe_parameters());

    chargefw::test::assert_calculation_provenance(charge_set, "mpeoe", "MPEOE_original");
    assert(charge_set.size() == 1);

    const auto& charges = charge_set.assignment(0).charges;

    assert(charges.size() == 5);

    const auto carbon_charge = charges[0];
    const auto hydrogen_charge = charges[1];

    assert(carbon_charge < 0.0);
    assert(hydrogen_charge > 0.0);

    assert(std::abs(charges[1] - charges[2]) < 1.0e-12);
    assert(std::abs(charges[1] - charges[3]) < 1.0e-12);
    assert(std::abs(charges[1] - charges[4]) < 1.0e-12);

    assert(std::abs(charges.total()) < 1.0e-12);

    return 0;
}
