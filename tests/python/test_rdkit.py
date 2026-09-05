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


class FakeBond:
    def GetBeginAtomIdx(self) -> int:
        return 0

    def GetEndAtomIdx(self) -> int:
        return 1

    def GetBondTypeAsDouble(self) -> float:
        return 1.0


class FakeMol:
    def __init__(self) -> None:
        self.atoms = (FakeAtom(0, 8, "O"), FakeAtom(1, 1, "H"))
        self.conformers = (FakeConformer(),)
        self.properties: dict[str, str] = {}

    def GetAtoms(self) -> tuple[FakeAtom, ...]:
        return self.atoms

    def GetBonds(self) -> tuple[FakeBond, ...]:
        return (FakeBond(),)

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

    @staticmethod
    def CreateAtomDoublePropertyList(molecule: FakeMol, name: str) -> None:
        molecule.properties[f"atom.dprop.{name}"] = "created"


class RdkitAdapterTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
