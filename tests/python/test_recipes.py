"""Executable documentation recipe checks."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import gemmi
import pytest

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


WATERS_SDF = f"{WATER_MOL}$$$$\n{WATER_MOL}$$$$\n"


def run_recipe(name: str, *arguments: str | Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(RECIPES / name), *map(str, arguments)],
        check=True,
        capture_output=True,
        text=True,
    )


def test_calculate_file_recipe(tmp_path: Path) -> None:
    input_path = tmp_path / "waters.sdf"
    input_path.write_text(WATERS_SDF, encoding="utf-8")

    completed = run_recipe("calculate_file.py", input_path, "--format", "sdf", "--method", "eem")

    lines = completed.stdout.splitlines()
    assert len(lines) == 2
    assert lines[0].startswith("molecule=0 conformer=0: [")
    assert lines[1].startswith("molecule=1 conformer=0: [")


def test_inspect_molecules_recipe(tmp_path: Path) -> None:
    input_path = tmp_path / "water.sdf"
    input_path.write_text(f"{WATER_MOL}$$$$\n", encoding="utf-8")

    run_recipe("inspect_molecules.py", input_path, "--format", "sdf")


def test_calculate_sdf_collection_recipe(tmp_path: Path) -> None:
    input_path = tmp_path / "waters.sdf"
    output_path = tmp_path / "result.json"
    input_path.write_text(WATERS_SDF, encoding="utf-8")

    run_recipe("calculate_sdf_collection.py", input_path, output_path)

    result = json.loads(output_path.read_text(encoding="utf-8"))
    assert result["status"] == "success"
    assert len(result["results"]) == 2


def test_gemmi_document_recipe_preserves_unchanged_mmcif(tmp_path: Path) -> None:
    input_path = tmp_path / "water.cif"
    output_path = tmp_path / "charged.cif"
    document = gemmi.read_pdb_string(WATER_PDB).make_mmcif_document()
    document.sole_block().set_pair("_audit.creation_method", "recipe-test")
    document.write_file(str(input_path))

    run_recipe("charge_gemmi_document.py", input_path, output_path, "--format", "mmcif")

    block = gemmi.cif.read_file(str(output_path)).sole_block()
    assert block.find_value("_audit.creation_method") == "recipe-test"
    assert "_sb_ncbr_partial_atomic_charges." in block.get_mmcif_category_names()


def test_fixed_ion_recipe(tmp_path: Path) -> None:
    input_path = PROJECT_ROOT / "tests" / "fixtures" / "synthetic" / "cif" / "fixed_ions.cif"
    output_path = tmp_path / "charges.json"

    run_recipe(
        "calculate_with_fixed_ions.py",
        input_path,
        "--format",
        "mmcif",
        "--ion",
        "MG",
        "--ion",
        "CA",
        "--result-json",
        output_path,
    )

    result = json.loads(output_path.read_text(encoding="utf-8"))
    assert result["status"] == "success"
    assert result["calculation_provenance"]["effective"]["fixed_ions"][0]["component_id"] == "MG"


def test_rdkit_conformer_recipe() -> None:
    pytest.importorskip("rdkit")
    run_recipe("analyze_rdkit_conformer_charges.py", "CCO", "--conformers", "3")


def test_compare_parameter_sets_recipe(tmp_path: Path) -> None:
    input_path = tmp_path / "water.sdf"
    input_path.write_text(f"{WATER_MOL}$$$$\n", encoding="utf-8")

    run_recipe("compare_parameter_sets.py", input_path, "--format", "sdf", "--method", "qeq")
