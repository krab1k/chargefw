"""Owned molecule and collection model checks."""

import unittest
from typing import Any, cast

import chargefw
import numpy as np


class MoleculeTests(unittest.TestCase):
    def test_molecule_owns_normalized_input(self) -> None:
        atomic_source = np.array([8, 1, 1, 6], dtype=np.int8)
        formal = np.array([0, 0, 0, -1], dtype=np.int8)
        bonds = np.array([[0, 1, 1], [0, 2, 1], [0, 3, 2]], dtype=np.int32)
        coordinates_source = np.arange(24, dtype=np.float32).reshape(2, 4, 3)
        atom_ids = ["O", 17, ("C", 1), "H"]
        conformer_ids = ["model-a", "model-b"]

        molecule = chargefw.Molecule(
            atomic_source,
            formal_charges=formal,
            bonds=bonds,
            coordinates=coordinates_source[:, :, ::-1],
            name="water-like",
            atom_names=["O", "H1", "H2", "C"],
            conformer_names=["a", "b"],
            source_name="fixture.sdf",
            record_index=3,
            record_id="record-4",
            atom_ids=atom_ids,
            conformer_ids=conformer_ids,
        )

        atomic_source[:] = 1
        formal[:] = 9
        bonds[:] = 0
        coordinates_source[:] = 99

        np.testing.assert_array_equal(molecule.atomic_numbers, [8, 1, 1, 6])
        np.testing.assert_array_equal(molecule.formal_charges, [0, 0, 0, -1])
        np.testing.assert_array_equal(molecule.bonds, [[0, 1, 1], [0, 2, 1], [0, 3, 2]])
        self.assertEqual(molecule.coordinates.shape, (2, 4, 3))
        self.assertEqual(molecule.coordinates[0, 0].tolist(), [2.0, 1.0, 0.0])
        self.assertEqual(molecule.source, chargefw.SourceIdentity("fixture.sdf", 3, "record-4"))
        self.assertEqual(molecule.atom_ids, tuple(atom_ids))
        self.assertEqual(molecule.conformer_ids, tuple(conformer_ids))
        self.assertEqual(molecule._native.atom_count, 4)
        self.assertEqual(molecule._native.bond_count, 3)
        self.assertEqual(molecule._native.conformer_count, 2)
        self.assertEqual(molecule._native.name, "water-like")

        readonly = molecule.atomic_numbers
        self.assertIs(readonly, molecule.atomic_numbers)
        with self.assertRaises(ValueError):
            readonly[0] = 2
        with self.assertRaises(ValueError):
            readonly.setflags(write=True)
        with self.assertRaises(ValueError):
            molecule._formal_charges[0] = 2
        with self.assertRaises(ValueError):
            molecule._bonds[0, 0] = 1
        with self.assertRaises(ValueError):
            molecule._coordinates[0, 0, 0] = 1.0
        with self.assertRaises(AttributeError):
            setattr(molecule, "name", "changed")
        self.assertEqual(
            repr(molecule),
            "Molecule(atom_count=4, bond_count=3, conformer_count=2, name='water-like')",
        )

    def test_collection_is_an_immutable_sequence(self) -> None:
        molecule = chargefw.Molecule([1])
        collection = chargefw.MoleculeCollection([molecule], name="fixture")

        self.assertEqual(len(collection), 1)
        self.assertIs(collection[0], molecule)
        self.assertEqual(collection[:], (molecule,))
        self.assertEqual(tuple(collection), (molecule,))
        self.assertEqual(collection.molecules, (molecule,))
        self.assertEqual(collection.name, "fixture")
        self.assertEqual(collection._native_molecules, (molecule._native,))
        self.assertEqual(repr(collection), "MoleculeCollection(molecules=1, name='fixture')")

    def test_coordinate_defaults_and_empty_molecule(self) -> None:
        no_coordinates = chargefw.Molecule([1])
        self.assertEqual(no_coordinates.coordinates.shape, (0, 1, 3))
        self.assertEqual(no_coordinates.conformer_count, 0)
        self.assertFalse(no_coordinates.has_coordinates)
        self.assertEqual(chargefw.Molecule([]).atomic_numbers.dtype, np.dtype(np.int64))

        one_conformer = chargefw.Molecule([1], coordinates=[[0.0, 0.0, 0.0]])
        self.assertEqual(one_conformer.coordinates.shape, (1, 1, 3))
        self.assertEqual(one_conformer.conformer_count, 1)

    def test_integer_iterables_and_empty_bonds_are_normalized(self) -> None:
        molecule = chargefw.Molecule(
            (atomic_number for atomic_number in [8, 1, 1]), bonds=[]
        )
        np.testing.assert_array_equal(molecule.atomic_numbers, [8, 1, 1])
        self.assertEqual(molecule.bonds.shape, (0, 3))

    def test_object_integer_arrays_are_range_checked(self) -> None:
        molecule = chargefw.Molecule(np.array([1, 8], dtype=object))
        np.testing.assert_array_equal(molecule.atomic_numbers, [1, 8])
        with self.assertRaises(ValueError):
            chargefw.Molecule(np.array([2**100], dtype=object))

    def test_invalid_inputs_are_rejected(self) -> None:
        invalid_cases = (
            (TypeError, lambda: chargefw.Molecule([1.0])),
            (ValueError, lambda: chargefw.Molecule([0])),
            (ValueError, lambda: chargefw.Molecule([1, 1], bonds=[[0, 1, 4]])),
            (
                ValueError,
                lambda: chargefw.Molecule([1, 1], bonds=[[0, 1, 1], [1, 0, 1]]),
            ),
            (ValueError, lambda: chargefw.Molecule([1], coordinates=[[np.nan, 0, 0]])),
            (
                ValueError,
                lambda: chargefw.Molecule([1], coordinates=np.empty((0, 1, 2))),
            ),
            (ValueError, lambda: chargefw.Molecule([1], atom_ids=cast(Any, [[]]))),
            (TypeError, lambda: chargefw.Molecule([1], atom_names="H")),
            (ValueError, lambda: chargefw.SourceIdentity(record_index=-1)),
            (TypeError, lambda: chargefw.SourceIdentity(record_index=cast(Any, 1.5))),
            (TypeError, lambda: chargefw.SourceIdentity(record_index=True)),
            (
                TypeError,
                lambda: chargefw.Molecule(
                    [1], record_index=cast(Any, np.bool_(True))
                ),
            ),
        )
        for error_type, operation in invalid_cases:
            with self.subTest(error_type=error_type, operation=operation), self.assertRaises(
                error_type
            ):
                operation()


if __name__ == "__main__":
    unittest.main()
