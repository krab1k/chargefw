#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/spatial_fragment.h>
#include <chargefw/parameters/classification/parameter_classification.h>

#include <cassert>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace parameters = chargefw::parameters;

namespace {

auto make_molecule() -> core::Molecule {
    auto atoms = std::vector{core::Atom{6, 0, "C0"}, core::Atom{8, -1, "O1"},
                             core::Atom{7, 1, "N2"}, core::Atom{1, 0, "H3"}};
    auto bonds = std::vector{core::Bond{0, 1, core::BondOrder::DOUBLE},
                             core::Bond{1, 2, core::BondOrder::SINGLE},
                             core::Bond{2, 3, core::BondOrder::SINGLE}};
    auto conformers =
        std::vector{core::Conformer{{core::Position{.x = 0.0}, core::Position{.x = 1.0},
                                     core::Position{.x = 2.0}, core::Position{.x = 10.0}},
                                    "model-1"},
                    core::Conformer{{core::Position{.x = 0.0}, core::Position{.x = 2.0},
                                     core::Position{.x = 4.0}, core::Position{.x = 20.0}},
                                    "model-2"}};
    return core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers), "fragment"};
}

template <typename Exception, typename Callable> auto throws(Callable&& callable) -> bool {
    try {
        std::forward<Callable>(callable)();
    } catch (const Exception&) {
        return true;
    }

    return false;
}

} // namespace

auto main() -> int {
    const auto molecule = make_molecule();
    const features::PreparedMolecule prepared{molecule};
    const auto fragment = features::build_spatial_fragment(prepared, 0, 1, 1.0);

    assert(fragment.molecule().name() == "fragment");
    assert(fragment.molecule().atom_count() == 3);
    assert(fragment.molecule().bond_count() == 2);
    assert(fragment.molecule().conformer_count() == 1);
    assert(fragment.molecule().conformer(0).name() == "model-1");
    assert(fragment.source_conformer_index() == 0);
    assert(fragment.center_local_atom_index() == 1);

    assert(std::ranges::equal(fragment.local_to_source_atom_indices(),
                              std::vector<std::size_t>{0, 1, 2}));
    assert(std::ranges::equal(fragment.source_to_local_atom_indices(),
                              std::vector<std::size_t>{0, 1, 2, features::no_source_index}));
    assert(std::ranges::equal(fragment.local_to_source_bond_indices(),
                              std::vector<std::size_t>{0, 1}));

    assert(fragment.molecule().atom(1).atomic_number() == 8);
    assert(fragment.molecule().atom(1).formal_charge() == -1);
    assert(fragment.molecule().atom(1).name() == "O1");
    assert(fragment.molecule().bond(0).first_atom_index() == 0);
    assert(fragment.molecule().bond(0).second_atom_index() == 1);
    assert(fragment.molecule().bond(0).order() == core::BondOrder::DOUBLE);
    assert(fragment.molecule().bond(1).first_atom_index() == 1);
    assert(fragment.molecule().bond(1).second_atom_index() == 2);
    assert(fragment.molecule().conformer(0)[0].x == 0.0);
    assert(fragment.molecule().conformer(0)[1].x == 1.0);
    assert(fragment.molecule().conformer(0)[2].x == 2.0);

    const parameters::ParameterClassification classification{
        parameters::AtomParameterClassification{{10, 11, 12, 13}},
        parameters::BondParameterClassification{{20, 21, 22}}};
    const auto projected = features::project_classification(classification, fragment);
    assert(std::ranges::equal(projected.atom().parameter_entry_indices(),
                              std::vector<std::size_t>{10, 11, 12}));
    assert(std::ranges::equal(projected.bond().parameter_entry_indices(),
                              std::vector<std::size_t>{20, 21}));

    const auto atom_only_projected = features::project_classification(
        parameters::ParameterClassification{
            parameters::AtomParameterClassification{{10, 11, 12, 13}}},
        fragment);
    assert(std::ranges::equal(atom_only_projected.atom().parameter_entry_indices(),
                              std::vector<std::size_t>{10, 11, 12}));
    assert(atom_only_projected.bond().empty());

    const auto bond_only_projected = features::project_classification(
        parameters::ParameterClassification{parameters::AtomParameterClassification{},
                                            parameters::BondParameterClassification{{20, 21, 22}}},
        fragment);
    assert(bond_only_projected.atom().empty());
    assert(std::ranges::equal(bond_only_projected.bond().parameter_entry_indices(),
                              std::vector<std::size_t>{20, 21}));

    const auto boundary_fragment = features::build_spatial_fragment(prepared, 0, 2, 8.0);
    assert(boundary_fragment.molecule().atom_count() == 4);
    assert(boundary_fragment.center_local_atom_index() == 2);
    assert(std::ranges::equal(boundary_fragment.local_to_source_atom_indices(),
                              std::vector<std::size_t>{0, 1, 2, 3}));

    const auto second_conformer_fragment = features::build_spatial_fragment(prepared, 1, 1, 2.0);
    assert(second_conformer_fragment.molecule().conformer(0).name() == "model-2");
    assert(second_conformer_fragment.molecule().conformer(0)[2].x == 4.0);

    const auto repeated_fragment = features::build_spatial_fragment(prepared, 0, 1, 1.0);
    assert(std::ranges::equal(repeated_fragment.local_to_source_atom_indices(),
                              fragment.local_to_source_atom_indices()));

    const core::Molecule single_atom_molecule{
        std::vector{core::Atom{1, 0, "H"}},
        {},
        std::vector{core::Conformer{{core::Position{.x = 3.0}}, "single"}},
        "single"};
    const features::PreparedMolecule single_atom_prepared{single_atom_molecule};
    const auto single_atom_fragment =
        features::build_spatial_fragment(single_atom_prepared, 0, 0, 1.0);
    assert(single_atom_fragment.molecule().atom_count() == 1);
    assert(single_atom_fragment.molecule().bond_count() == 0);
    assert(single_atom_fragment.center_local_atom_index() == 0);

    const core::Molecule empty_molecule{{}, {}, {core::Conformer{{}, "empty"}}, "empty"};
    const features::PreparedMolecule empty_prepared{empty_molecule};
    assert(throws<std::out_of_range>([&empty_prepared] {
        static_cast<void>(features::build_spatial_fragment(empty_prepared, 0, 0, 1.0));
    }));
    assert(throws<std::out_of_range>(
        [&prepared] { static_cast<void>(features::build_spatial_fragment(prepared, 2, 1, 1.0)); }));
    assert(throws<std::out_of_range>(
        [&prepared] { static_cast<void>(features::build_spatial_fragment(prepared, 0, 4, 1.0)); }));
    assert(throws<std::invalid_argument>(
        [&prepared] { static_cast<void>(features::build_spatial_fragment(prepared, 0, 1, 0.0)); }));
    assert(throws<std::invalid_argument>([&prepared] {
        static_cast<void>(features::build_spatial_fragment(prepared, 0, 1, -1.0));
    }));
    assert(throws<std::invalid_argument>([&prepared] {
        static_cast<void>(features::build_spatial_fragment(
            prepared, 0, 1, std::numeric_limits<double>::quiet_NaN()));
    }));
    assert(throws<std::invalid_argument>([&prepared] {
        static_cast<void>(features::build_spatial_fragment(
            prepared, 0, 1, std::numeric_limits<double>::infinity()));
    }));

    return 0;
}
