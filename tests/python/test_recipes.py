"""Executable documentation recipe checks."""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import gemmi

PROJECT_ROOT = Path(__file__).parents[2]
RECIPES = PROJECT_ROOT / "docs" / "recipes"

WATER_PDB = """\
HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O
HETATM    2  H1  HOH A   1       0.957   0.000   0.000  1.00 20.00           H
HETATM    3  H2  HOH A   1      -0.240   0.927   0.000  1.00 20.00           H
CONECT    1    2    3
END
"""

WATER_MOL = """water
  ChargeFW

  3  2  0  0  0  0  0  0  0  0  1 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9570    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2400    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
"""


class RecipeTests(unittest.TestCase):
    def test_calculate_file_recipe(self) -> None:
        with TemporaryDirectory() as directory:
            input_path = Path(directory) / "water.sdf"
            input_path.write_text(f"{WATER_MOL}$$$$\n", encoding="utf-8")

            completed = subprocess.run(
                [
                    sys.executable,
                    str(RECIPES / "calculate_file.py"),
                    str(input_path),
                    "--format",
                    "sdf",
                    "--method",
                    "qeq",
                    "--parameter-set",
                    "QEq_original",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.stdout.count("["), 1)

    def test_inspect_molecules_recipe(self) -> None:
        with TemporaryDirectory() as directory:
            input_path = Path(directory) / "water.sdf"
            input_path.write_text(f"{WATER_MOL}$$$$\n", encoding="utf-8")

            completed = subprocess.run(
                [
                    sys.executable,
                    str(RECIPES / "inspect_molecules.py"),
                    str(input_path),
                    "--format",
                    "sdf",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertIn("1 molecule(s)", completed.stdout)
            self.assertIn("atoms=3 bonds=2 conformers=1 formal_charge=0", completed.stdout)

    def test_calculate_sdf_collection_recipe(self) -> None:
        with TemporaryDirectory() as directory:
            input_path = Path(directory) / "waters.sdf"
            output_path = Path(directory) / "result.json"
            input_path.write_text(f"{WATER_MOL}$$$$\n{WATER_MOL}$$$$\n", encoding="utf-8")

            subprocess.run(
                [
                    sys.executable,
                    str(RECIPES / "calculate_sdf_collection.py"),
                    str(input_path),
                    str(output_path),
                ],
                check=True,
            )
            result = json.loads(output_path.read_text(encoding="utf-8"))
            self.assertEqual(result["status"], "success")
            self.assertEqual(len(result["results"]), 2)

    def test_gemmi_document_recipe_preserves_mmcif(self) -> None:
        with TemporaryDirectory() as directory:
            input_path = Path(directory) / "water.cif"
            output_path = Path(directory) / "charged.cif"
            structure = gemmi.read_pdb_string(WATER_PDB)
            document = structure.make_mmcif_document()
            document.sole_block().set_pair("_audit.creation_method", "recipe-test")
            document.write_file(str(input_path))

            subprocess.run(
                [
                    sys.executable,
                    str(RECIPES / "charge_gemmi_document.py"),
                    str(input_path),
                    str(output_path),
                    "--format",
                    "mmcif",
                ],
                check=True,
            )

            document = gemmi.cif.read_file(str(output_path))
            self.assertEqual(
                document.sole_block().find_value("_audit.creation_method"), "recipe-test"
            )
            self.assertIn(
                "_sb_ncbr_partial_atomic_charges.", document.sole_block().get_mmcif_category_names()
            )

    @unittest.skipUnless(importlib.util.find_spec("rdkit"), "RDKit is not installed")
    def test_rdkit_conformer_recipe(self) -> None:
        completed = subprocess.run(
            [
                sys.executable,
                str(RECIPES / "analyze_rdkit_conformer_charges.py"),
                "CCO",
                "--conformers",
                "3",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertIn("Method: qeq/QEq_original", completed.stdout)
        self.assertIn("equally weighted, unoptimized ETKDGv3 structures", completed.stdout)
        self.assertIn("Atom          mean         std", completed.stdout)

    def test_compare_parameter_sets_recipe(self) -> None:
        with TemporaryDirectory() as directory:
            input_path = Path(directory) / "water.sdf"
            input_path.write_text(f"{WATER_MOL}$$$$\n", encoding="utf-8")

            completed = subprocess.run(
                [
                    sys.executable,
                    str(RECIPES / "compare_parameter_sets.py"),
                    str(input_path),
                    "--format",
                    "sdf",
                    "--method",
                    "qeq",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertIn("Reference: QEq_original", completed.stdout)
            self.assertIn("RMS difference", completed.stdout)


if __name__ == "__main__":
    unittest.main()
