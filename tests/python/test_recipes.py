"""Executable documentation recipe checks."""

from __future__ import annotations

import importlib.util
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


class RecipeTests(unittest.TestCase):
    def test_gemmi_structure_recipe(self) -> None:
        with TemporaryDirectory() as directory:
            input_path = Path(directory) / "water.pdb"
            output_path = Path(directory) / "charged.cif"
            input_path.write_text(WATER_PDB, encoding="utf-8")

            subprocess.run(
                [
                    sys.executable,
                    str(RECIPES / "charge_structure_with_gemmi.py"),
                    str(input_path),
                    str(output_path),
                    "--method",
                    "formal",
                ],
                check=True,
            )

            document = gemmi.cif.read_file(str(output_path))
            self.assertIn(
                "_sb_ncbr_partial_atomic_charges.", document.sole_block().get_mmcif_category_names()
            )

    @unittest.skipUnless(importlib.util.find_spec("rdkit"), "RDKit is not installed")
    def test_rdkit_conformer_recipe(self) -> None:
        completed = subprocess.run(
            [
                sys.executable,
                str(RECIPES / "conformational_charge_hotspots.py"),
                "CCO",
                "--conformers",
                "3",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertIn("Lowest charge variance", completed.stdout)
        self.assertIn("Highest charge variance", completed.stdout)


if __name__ == "__main__":
    unittest.main()
