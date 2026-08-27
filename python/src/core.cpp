#include "bindings.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>

#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;

namespace chargefw::python {
namespace {

using integer_array_1d =
    nb::ndarray<const std::int64_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using integer_array_2d =
    nb::ndarray<const std::int64_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using coordinate_array = nb::ndarray<const double, nb::ndim<3>, nb::c_contig, nb::device::cpu>;

auto sequence_string(const nb::sequence& values, const std::size_t index) -> std::string {
    return nb::cast<std::string>(values[index]);
}

auto make_molecule(integer_array_1d atomic_numbers, integer_array_1d formal_charges,
                   integer_array_2d bonds, coordinate_array coordinates, nb::sequence atom_names,
                   nb::sequence conformer_names, std::string name) -> core::Molecule {
    const auto atom_count = static_cast<std::size_t>(atomic_numbers.shape(0));

    if (formal_charges.shape(0) != atomic_numbers.shape(0)) {
        throw std::invalid_argument{"formal charge count must match atomic number count"};
    }
    if (bonds.shape(1) != 3) {
        throw std::invalid_argument{"bonds must have three columns"};
    }
    if (coordinates.shape(1) != static_cast<std::int64_t>(atom_count) ||
        coordinates.shape(2) != 3) {
        throw std::invalid_argument{"coordinate shape does not match the molecule"};
    }

    std::vector<core::Atom> atoms;
    atoms.reserve(atom_count);
    for (std::size_t index = 0; index < atom_count; ++index) {
        atoms.emplace_back(static_cast<int>(atomic_numbers.data()[index]),
                           static_cast<int>(formal_charges.data()[index]),
                           sequence_string(atom_names, index));
    }

    const auto bond_count = static_cast<std::size_t>(bonds.shape(0));
    std::vector<core::Bond> native_bonds;
    native_bonds.reserve(bond_count);
    for (std::size_t index = 0; index < bond_count; ++index) {
        const auto offset = index * 3;
        native_bonds.emplace_back(
            static_cast<std::size_t>(bonds.data()[offset]),
            static_cast<std::size_t>(bonds.data()[offset + 1]),
            core::bond_order_from_value(static_cast<int>(bonds.data()[offset + 2])));
    }

    const auto conformer_count = static_cast<std::size_t>(coordinates.shape(0));
    std::vector<core::Conformer> conformers;
    conformers.reserve(conformer_count);
    for (std::size_t conformer_index = 0; conformer_index < conformer_count; ++conformer_index) {
        std::vector<core::Position> positions;
        positions.reserve(atom_count);
        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            const auto offset = (conformer_index * atom_count + atom_index) * 3;
            positions.push_back(core::Position{
                .x = coordinates.data()[offset],
                .y = coordinates.data()[offset + 1],
                .z = coordinates.data()[offset + 2],
            });
        }
        conformers.emplace_back(std::move(positions),
                                sequence_string(conformer_names, conformer_index));
    }

    return core::Molecule{std::move(atoms), std::move(native_bonds), std::move(conformers),
                          std::move(name)};
}

auto make_collection(nb::sequence molecules, std::string name) -> core::MoleculeCollection {
    const auto molecule_count = nb::len(molecules);
    std::vector<core::Molecule> native_molecules;
    native_molecules.reserve(molecule_count);
    for (std::size_t index = 0; index < molecule_count; ++index) {
        native_molecules.emplace_back(nb::cast<const core::Molecule&>(molecules[index]));
    }
    return core::MoleculeCollection{std::move(native_molecules), std::move(name)};
}

} // namespace

void bind_core(nb::module_& module) {
    nb::class_<core::Molecule>(module, "_NativeMolecule")
        .def_prop_ro("name", &core::Molecule::name)
        .def_prop_ro("atom_count", &core::Molecule::atom_count)
        .def_prop_ro("bond_count", &core::Molecule::bond_count)
        .def_prop_ro("conformer_count", &core::Molecule::conformer_count);
    nb::class_<core::MoleculeCollection>(module, "_NativeMoleculeCollection")
        .def_prop_ro("name", &core::MoleculeCollection::name)
        .def_prop_ro("size", &core::MoleculeCollection::size);

    module.def("_make_molecule", &make_molecule, nb::arg("atomic_numbers"),
               nb::arg("formal_charges"), nb::arg("bonds"), nb::arg("coordinates"),
               nb::arg("atom_names"), nb::arg("conformer_names"), nb::arg("name"));
    module.def("_make_collection", &make_collection, nb::arg("molecules"), nb::arg("name"));
}

} // namespace chargefw::python
