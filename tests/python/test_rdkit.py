"""Optional RDKit conversion, attachment, and serialization checks."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import chargefw
import numpy as np
from chargefw.io import rdkit as chargefw_rdkit

try:
    from rdkit import Chem  # type: ignore[import-not-found]
except ModuleNotFoundError:
    Chem = None


MOL_TEXT = """water
  ChargeFW

  2  1  0  0  0  0  0  0  0  0  1 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9600    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
"""


class FakePoint:
    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z


class FakeConformer:
    def GetId(self) -> int:
        return 0

    def GetAtomPosition(self, index: int) -> FakePoint:
        return (FakePoint(0.0, 0.0, 0.0), FakePoint(0.96, 0.0, 0.0))[index]


class FakeAtom:
    def __init__(self, index: int, atomic_number: int, symbol: str) -> None:
        self.index = index
        self.atomic_number = atomic_number
        self.symbol = symbol
        self.properties: dict[str, float] = {}

    def GetIdx(self) -> int:
        return self.index

    def GetAtomicNum(self) -> int:
        return self.atomic_number

    def GetFormalCharge(self) -> int:
        return 0

    def GetSymbol(self) -> str:
        return self.symbol

    def HasProp(self, name: str) -> bool:
        return name in self.properties

    def GetProp(self, name: str) -> str:
        raise KeyError(name)

    def SetDoubleProp(self, name: str, value: float) -> None:
        self.properties[name] = value


class FakeBondType:
    SINGLE = "SINGLE"
    DOUBLE = "DOUBLE"
    TRIPLE = "TRIPLE"
    AROMATIC = "AROMATIC"
    DATIVEONE = "DATIVEONE"
    DATIVE = "DATIVE"
    DATIVEL = "DATIVEL"
    DATIVER = "DATIVER"
    ZERO = "ZERO"


class FakeBond:
    def __init__(
        self, bond_type: str = FakeBondType.SINGLE, *, aromatic: bool = False, query: bool = False
    ) -> None:
        self.bond_type = bond_type
        self.aromatic = aromatic
        self.query = query

    def GetBeginAtomIdx(self) -> int:
        return 0

    def GetEndAtomIdx(self) -> int:
        return 1

    def HasQuery(self) -> bool:
        return self.query

    def GetBondType(self) -> str:
        return self.bond_type

    def GetIsAromatic(self) -> bool:
        return self.aromatic


class FakeMol:
    def __init__(self) -> None:
        self.atoms = (FakeAtom(0, 8, "O"), FakeAtom(1, 1, "H"))
        self.bonds = (FakeBond(),)
        self.conformers: tuple[FakeConformer, ...] = (FakeConformer(),)
        self.properties: dict[str, str] = {}

    def GetAtoms(self) -> tuple[FakeAtom, ...]:
        return self.atoms

    def GetBonds(self) -> tuple[FakeBond, ...]:
        return self.bonds

    def GetConformers(self) -> tuple[FakeConformer, ...]:
        return self.conformers

    def HasProp(self, name: str) -> bool:
        return name in self.properties

    def GetProp(self, name: str) -> str:
        return self.properties[name]

    def GetNumAtoms(self) -> int:
        return len(self.atoms)

    def GetAtomWithIdx(self, index: int) -> FakeAtom:
        return self.atoms[index]


class FakeChemistry:
    Mol = FakeMol
    BondType = FakeBondType

    @staticmethod
    def CreateAtomDoublePropertyList(molecule: FakeMol, name: str) -> None:
        molecule.properties[f"atom.dprop.{name}"] = "created"


class RdkitAdapterTests(unittest.TestCase):
    def test_conversion_without_conformers_supports_coordinate_independent_methods(self) -> None:
        target = FakeMol()
        target.conformers = ()

        with patch.object(chargefw_rdkit, "_require_rdkit", return_value=FakeChemistry):
            molecule = chargefw_rdkit.from_mol(target)

        self.assertEqual(molecule.coordinates.shape, (0, 2, 3))
        self.assertFalse(molecule.has_coordinates)
        result = chargefw.calculate(molecule, method="formal")
        np.testing.assert_array_equal(result.assignments[0].values, [0.0, 0.0])

    def test_conversion_and_charge_attachment(self) -> None:
        target = FakeMol()
        with patch.object(chargefw_rdkit, "_require_rdkit", return_value=FakeChemistry):
            molecule = chargefw_rdkit.from_mol(target, source_name="water")
            result = chargefw.calculate(molecule, method="formal")
            chargefw_rdkit.attach_charges(target, result)
            with self.assertRaisesRegex(ValueError, "already has property"):
                chargefw_rdkit.attach_charges(target, result)
            chargefw_rdkit.attach_charges(target, result, overwrite=True)

        self.assertEqual(molecule.atom_ids, (0, 1))
        self.assertEqual(target.atoms[0].properties["ChargeFWPartialCharge"], 0.0)
        self.assertEqual(target.properties["atom.dprop.ChargeFWPartialCharge"], "created")

    def test_conversion_rejects_unsupported_bond_types_by_default(self) -> None:
        for bond_type in (FakeBondType.AROMATIC, FakeBondType.DATIVE, FakeBondType.DATIVEONE):
            with self.subTest(bond_type=bond_type):
                target = FakeMol()
                target.bonds = (FakeBond(bond_type),)
                with (
                    patch.object(chargefw_rdkit, "_require_rdkit", return_value=FakeChemistry),
                    self.assertRaisesRegex(ValueError, f"unsupported bond type {bond_type}"),
                ):
                    chargefw_rdkit.from_mol(target)

    def test_single_conversion_normalizes_only_aromatic_and_dative_bonds(self) -> None:
        for bond_type in (
            FakeBondType.AROMATIC,
            FakeBondType.DATIVEONE,
            FakeBondType.DATIVE,
            FakeBondType.DATIVEL,
            FakeBondType.DATIVER,
        ):
            with self.subTest(bond_type=bond_type):
                target = FakeMol()
                target.bonds = (FakeBond(bond_type),)
                with patch.object(chargefw_rdkit, "_require_rdkit", return_value=FakeChemistry):
                    molecule = chargefw_rdkit.from_mol(target, bond_conversion="single")
                self.assertEqual(molecule.bonds.tolist(), [[0, 1, 1]])

        target = FakeMol()
        target.bonds = (FakeBond(FakeBondType.ZERO),)
        with (
            patch.object(chargefw_rdkit, "_require_rdkit", return_value=FakeChemistry),
            self.assertRaisesRegex(ValueError, "unsupported bond type ZERO"),
        ):
            chargefw_rdkit.from_mol(target, bond_conversion="single")

        target.bonds = (FakeBond(query=True),)
        with (
            patch.object(chargefw_rdkit, "_require_rdkit", return_value=FakeChemistry),
            self.assertRaisesRegex(ValueError, "query bonds cannot be imported"),
        ):
            chargefw_rdkit.from_mol(target, bond_conversion="single")

    def test_bond_conversion_is_validated(self) -> None:
        with (
            patch.object(chargefw_rdkit, "_require_rdkit", return_value=FakeChemistry),
            self.assertRaisesRegex(TypeError, "bond_conversion must be a string"),
        ):
            chargefw_rdkit.from_mol(FakeMol(), bond_conversion=True)  # type: ignore[arg-type]
        with (
            patch.object(chargefw_rdkit, "_require_rdkit", return_value=FakeChemistry),
            self.assertRaisesRegex(ValueError, "bond_conversion must be 'none' or 'single'"),
        ):
            chargefw_rdkit.from_mol(FakeMol(), bond_conversion="all")  # type: ignore[arg-type]

    def test_attachment_requires_a_bijective_atom_mapping(self) -> None:
        target = FakeMol()
        target.atoms = (FakeAtom(0, 8, "O"), FakeAtom(1, 8, "O"))
        molecule = chargefw.Molecule([8, 8], atom_ids=[0, 0])
        result = chargefw.calculate(molecule, method="formal")

        with (
            patch.object(chargefw_rdkit, "_require_rdkit", return_value=FakeChemistry),
            self.assertRaisesRegex(ValueError, "map each target atom exactly once"),
        ):
            chargefw_rdkit.attach_charges(target, result)

        self.assertFalse(any(atom.HasProp("ChargeFWPartialCharge") for atom in target.atoms))

    def test_attachment_accepts_index_compatible_atom_ids(self) -> None:
        target = FakeMol()
        molecule = chargefw.Molecule([1, 8], atom_ids=[np.int64(1), np.int64(0)])
        result = chargefw.calculate(molecule, method="formal")

        with patch.object(chargefw_rdkit, "_require_rdkit", return_value=FakeChemistry):
            chargefw_rdkit.attach_charges(target, result)

        self.assertTrue(all(atom.HasProp("ChargeFWPartialCharge") for atom in target.atoms))

    def test_missing_dependency_is_actionable(self) -> None:
        error = ModuleNotFoundError("No module named 'rdkit'")
        error.name = "rdkit"
        with (
            patch.object(chargefw_rdkit, "import_module", side_effect=error),
            self.assertRaisesRegex(ImportError, "independently installed RDKit"),
        ):
            chargefw_rdkit.from_mol(object())

    @unittest.skipIf(Chem is None, "RDKit is not installed")
    def test_real_rdkit_conversion_attachment_and_sd_serialization(self) -> None:
        assert Chem is not None
        target = Chem.MolFromMolBlock(MOL_TEXT, sanitize=False, removeHs=False)
        self.assertIsNotNone(target)
        molecule = chargefw_rdkit.from_mol(target, source_name="water.mol")
        result = chargefw.calculate(molecule, method="formal")

        chargefw_rdkit.attach_charges(target, result)

        self.assertEqual(target.GetAtomWithIdx(0).GetDoubleProp("ChargeFWPartialCharge"), 0.0)
        with TemporaryDirectory() as directory:
            path = Path(directory) / "charged.sdf"
            writer = Chem.SDWriter(str(path))
            writer.write(target)
            writer.close()
            loaded = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)[0]
            self.assertIsNotNone(loaded)
            self.assertEqual(loaded.GetAtomWithIdx(0).GetDoubleProp("ChargeFWPartialCharge"), 0.0)

    @unittest.skipIf(Chem is None, "RDKit is not installed")
    def test_real_rdkit_aromatic_and_dative_bond_conversion(self) -> None:
        assert Chem is not None
        aromatic = Chem.MolFromSmiles("c1ccccc1")
        self.assertIsNotNone(aromatic)
        with self.assertRaisesRegex(ValueError, "unsupported bond type AROMATIC"):
            chargefw_rdkit.from_mol(aromatic)
        normalized = chargefw_rdkit.from_mol(aromatic, bond_conversion="single")
        self.assertTrue(np.all(normalized.bonds[:, 2] == 1))

        for name in ("DATIVE", "DATIVEONE"):
            with self.subTest(bond_type=name):
                editable = Chem.RWMol()
                editable.AddAtom(Chem.Atom(7))
                editable.AddAtom(Chem.Atom(26))
                editable.AddBond(0, 1, getattr(Chem.BondType, name))
                target = editable.GetMol()
                self.assertEqual(target.GetBondWithIdx(0).GetBondTypeAsDouble(), 1.0)
                with self.assertRaisesRegex(ValueError, f"unsupported bond type {name}"):
                    chargefw_rdkit.from_mol(target)
                normalized = chargefw_rdkit.from_mol(target, bond_conversion="single")
                self.assertEqual(normalized.bonds.tolist(), [[0, 1, 1]])


if __name__ == "__main__":
    unittest.main()
