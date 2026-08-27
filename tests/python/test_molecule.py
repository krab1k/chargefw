"""Owned molecule and collection model checks."""

import numpy as np

import chargefw


def raises(error_type, function):
    try:
        function()
    except error_type:
        return
    raise AssertionError(f"expected {error_type.__name__}")


atomic_source = np.array([8, 1, 1, 6], dtype=np.int8)
atomic = atomic_source[::1]
formal = np.array([0, 0, 0, -1], dtype=np.int8)
bonds = np.array([[0, 1, 1], [0, 2, 1], [0, 3, 2]], dtype=np.int32)
coordinates_source = np.arange(24, dtype=np.float32).reshape(2, 4, 3)
coordinates = coordinates_source[:, :, ::-1]
atom_ids = ["O", 17, ("C", 1), "H"]
conformer_ids = ["model-a", "model-b"]

molecule = chargefw.Molecule(
    atomic,
    formal_charges=formal,
    bonds=bonds,
    coordinates=coordinates,
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
assert np.array_equal(molecule.atomic_numbers, [8, 1, 1, 6])
assert np.array_equal(molecule.formal_charges, [0, 0, 0, -1])
assert np.array_equal(molecule.bonds, bonds * 0 + [[0, 1, 1], [0, 2, 1], [0, 3, 2]])
assert molecule.coordinates.shape == (2, 4, 3)
assert molecule.coordinates[0, 0].tolist() == [2.0, 1.0, 0.0]
assert molecule.source == chargefw.SourceIdentity("fixture.sdf", 3, "record-4")
assert molecule.source_atom_ids == tuple(atom_ids)
assert molecule.source_conformer_ids == tuple(conformer_ids)
assert molecule._native.atom_count == 4
assert molecule._native.bond_count == 3
assert molecule._native.conformer_count == 2

readonly = molecule.atomic_numbers
raises(ValueError, lambda: readonly.__setitem__(0, 2))
raises(AttributeError, lambda: setattr(molecule, "name", "changed"))

collection = chargefw.MoleculeCollection([molecule], name="fixture")
assert len(collection) == 1
assert collection[0] is molecule
assert collection.molecules == (molecule,)
assert collection.name == "fixture"
assert collection._native.size == 1

no_coordinates = chargefw.Molecule([1])
assert no_coordinates.coordinates.shape == (0, 1, 3)
assert no_coordinates.conformer_count == 0
assert not no_coordinates.has_coordinates
assert chargefw.Molecule([]).atomic_numbers.dtype == np.int64

raises(TypeError, lambda: chargefw.Molecule([1.0]))
raises(ValueError, lambda: chargefw.Molecule([0]))
raises(ValueError, lambda: chargefw.Molecule([1, 1], bonds=[[0, 1, 4]]))
raises(ValueError, lambda: chargefw.Molecule([1, 1], bonds=[[0, 1, 1], [1, 0, 1]]))
raises(ValueError, lambda: chargefw.Molecule([1], coordinates=[[np.nan, 0, 0]]))
raises(ValueError, lambda: chargefw.Molecule([1], coordinates=np.empty((0, 1, 2))))
raises(ValueError, lambda: chargefw.Molecule([1], atom_ids=[[]]))
raises(TypeError, lambda: chargefw.Molecule([1], atom_names="H"))
raises(ValueError, lambda: chargefw.SourceIdentity(record_index=-1))
raises(TypeError, lambda: chargefw.SourceIdentity(record_index=1.5))
