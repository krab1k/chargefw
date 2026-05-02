#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_calculation.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/atom_parameters.h>
#include <chargefw/parameters/common_parameters.h>
#include <chargefw/parameters/parameter_key.h>
#include <chargefw/parameters/parameter_set.h>
#include <chargefw/parameters/parameter_set_metadata.h>

#include <cassert>
#include <cmath>
#include <string_view>
#include <utility>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

auto atom_key(const int atomic_number) -> parameters::AtomParameterKey {
    return parameters::AtomParameterKey{.atomic_number = atomic_number,
                                        .classification =
                                            parameters::AtomParameterClassificationKind::PLAIN,
                                        .type = "*"};
}

auto make_hf() -> core::Molecule {
    std::vector atoms{core::Atom{1, 0, "H"}, core::Atom{9, 0, "F"}};

    std::vector bonds{core::Bond{0, 1, core::BondOrder::SINGLE}};

    return core::Molecule{std::move(atoms), std::move(bonds), {}, "hf"};
}

auto make_collection() -> core::MoleculeCollection {
    std::vector molecules{make_hf()};

    return core::MoleculeCollection{std::move(molecules), "test"};
}

auto make_peoe_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "peoe-test-parameters", .method_id = "peoe", .name = "PEOE test parameters"},
        parameters::CommonParameters{{{.name = "dampH", .value = 20.02}}},
        parameters::AtomParameters{{{.key = atom_key(1),
                                     .parameters = {{.name = "A", .value = 7.17},
                                                    {.name = "B", .value = 6.24},
                                                    {.name = "C", .value = -0.56}}},
                                    {.key = atom_key(9),
                                     .parameters = {{.name = "A", .value = 12.06},
                                                    {.name = "B", .value = 13.85},
                                                    {.name = "C", .value = 3.98}}}}}};
}

} // namespace

auto main() -> int {
    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared{collection};

    const auto& registry = methods::method_registry();
    const auto* peoe = registry.find("peoe");

    assert(peoe != nullptr);
    assert(peoe->requires_parameters());

    const std::vector<const methods::Method*> candidate_methods{peoe};

    const std::vector parameter_sets{make_peoe_parameters()};

    const auto applicability =
        methods::find_applicable_methods(prepared, candidate_methods, parameter_sets);

    assert(applicability.applicable.size() == 1);
    assert(applicability.rejected.empty());

    const auto charge_set = methods::calculate_charges(applicability.applicable.front(), prepared);

    assert(charge_set.method_id() == std::string_view{"peoe"});
    assert(charge_set.parameter_set_id().has_value());
    assert(*charge_set.parameter_set_id() == std::string_view{"peoe-test-parameters"});
    assert(charge_set.assignment_count() == 1);

    const auto& charges = charge_set.assignment(0).charges;

    assert(charges.size() == 2);

    const auto h_charge = charges[0];
    const auto f_charge = charges[1];

    assert(h_charge > 0.0);
    assert(f_charge < 0.0);
    assert(std::abs(charges.total()) < 1.0e-12);

    return 0;
}