#include <chargefw/config.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;

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
                   nb::sequence conformer_names, std::string name) -> chargefw::core::Molecule {
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

    std::vector<chargefw::core::Atom> atoms;
    atoms.reserve(atom_count);
    for (std::size_t index = 0; index < atom_count; ++index) {
        atoms.emplace_back(static_cast<int>(atomic_numbers.data()[index]),
                           static_cast<int>(formal_charges.data()[index]),
                           sequence_string(atom_names, index));
    }

    const auto bond_count = static_cast<std::size_t>(bonds.shape(0));
    std::vector<chargefw::core::Bond> native_bonds;
    native_bonds.reserve(bond_count);
    for (std::size_t index = 0; index < bond_count; ++index) {
        const auto offset = index * 3;
        native_bonds.emplace_back(
            static_cast<std::size_t>(bonds.data()[offset]),
            static_cast<std::size_t>(bonds.data()[offset + 1]),
            chargefw::core::bond_order_from_value(static_cast<int>(bonds.data()[offset + 2])));
    }

    const auto conformer_count = static_cast<std::size_t>(coordinates.shape(0));
    std::vector<chargefw::core::Conformer> conformers;
    conformers.reserve(conformer_count);
    for (std::size_t conformer_index = 0; conformer_index < conformer_count; ++conformer_index) {
        std::vector<chargefw::core::Position> positions;
        positions.reserve(atom_count);
        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            const auto offset = (conformer_index * atom_count + atom_index) * 3;
            positions.push_back(chargefw::core::Position{
                .x = coordinates.data()[offset],
                .y = coordinates.data()[offset + 1],
                .z = coordinates.data()[offset + 2],
            });
        }
        conformers.emplace_back(std::move(positions),
                                sequence_string(conformer_names, conformer_index));
    }

    return chargefw::core::Molecule{std::move(atoms), std::move(native_bonds),
                                    std::move(conformers), std::move(name)};
}

auto make_collection(nb::sequence molecules, std::string name)
    -> chargefw::core::MoleculeCollection {
    std::vector<chargefw::core::Molecule> native_molecules;
    native_molecules.reserve(nb::len(molecules));
    for (std::size_t index = 0; index < nb::len(molecules); ++index) {
        native_molecules.emplace_back(nb::cast<const chargefw::core::Molecule&>(molecules[index]));
    }
    return chargefw::core::MoleculeCollection{std::move(native_molecules), std::move(name)};
}

} // namespace

NB_MODULE(_chargefw, module) {
    module.doc() = "Private native extension for ChargeFW.";
    module.def("version", [] { return CHARGEFW_VERSION_STRING; });

    nb::class_<chargefw::core::Molecule>(module, "_NativeMolecule")
        .def_prop_ro("name", &chargefw::core::Molecule::name)
        .def_prop_ro("atom_count", &chargefw::core::Molecule::atom_count)
        .def_prop_ro("bond_count", &chargefw::core::Molecule::bond_count)
        .def_prop_ro("conformer_count", &chargefw::core::Molecule::conformer_count);
    nb::class_<chargefw::core::MoleculeCollection>(module, "_NativeMoleculeCollection")
        .def_prop_ro("name", &chargefw::core::MoleculeCollection::name)
        .def_prop_ro("size", &chargefw::core::MoleculeCollection::size);

    module.def("_make_molecule", &make_molecule, nb::arg("atomic_numbers"),
               nb::arg("formal_charges"), nb::arg("bonds"), nb::arg("coordinates"),
               nb::arg("atom_names"), nb::arg("conformer_names"), nb::arg("name"));
    module.def("_make_collection", &make_collection, nb::arg("molecules"), nb::arg("name"));
}
