#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/parameter_set_io.h>

#include <cassert>
#include <cmath>
#include <filesystem>
#include <string_view>
#include <utility>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

#ifndef CHARGEFW_TEST_PARAMETER_DIR
#error "CHARGEFW_TEST_PARAMETER_DIR must be defined"
#endif

auto make_methane() -> core::Molecule {
    std::vector atoms{
        core::Atom{6, 0, "C"},
        core::Atom{1, 0, "H1"},
        core::Atom{1, 0, "H2"},
        core::Atom{1, 0, "H3"},
        core::Atom{1, 0, "H4"}
    };

    std::vector bonds{
        core::Bond{0, 1, core::BondOrder::SINGLE},
        core::Bond{0, 2, core::BondOrder::SINGLE},
        core::Bond{0, 3, core::BondOrder::SINGLE},
        core::Bond{0, 4, core::BondOrder::SINGLE}
    };

    return core::Molecule{
        std::move(atoms),
        std::move(bonds),
        {},
        "methane"
    };
}

auto make_collection() -> core::MoleculeCollection {
    std::vector molecules{
        make_methane()
    };

    return core::MoleculeCollection{
        std::move(molecules),
        "test"
    };
}

auto load_mpeoe_parameters() -> std::vector<parameters::ParameterSet> {
    const auto path =
        std::filesystem::path{CHARGEFW_TEST_PARAMETER_DIR} / "MPEOE_original.json";

    std::vector<parameters::ParameterSet> parameter_sets;
    parameter_sets.push_back(parameters::load_parameter_set_json_file(path));
    return parameter_sets;
}

} // namespace

auto main() -> int {
    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared{collection};

    const auto& registry = methods::method_registry();
    const auto* mpeoe = registry.find("mpeoe");

    assert(mpeoe != nullptr);
    assert(mpeoe->requires_parameters());

    const std::vector<const methods::Method*> candidate_methods{
        mpeoe
    };

    const auto parameter_sets = load_mpeoe_parameters();

    const auto applicability = methods::find_applicable_methods(
        prepared,
        candidate_methods,
        parameter_sets
    );

    assert(applicability.applicable.size() == 1);
    assert(applicability.rejected.empty());

    const auto charge_set = methods::calculate_charges(
        applicability.applicable.front(),
        prepared
    );

    assert(charge_set.method_id() == std::string_view{"mpeoe"});
    assert(charge_set.parameter_set_id().has_value());
    assert(*charge_set.parameter_set_id() == std::string_view{"MPEOE_original"});
    assert(charge_set.assignment_count() == 1);

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