"""Generated molecular and calculation-result output checks."""

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, cast

import chargefw
import gemmi


def water(*, conformers: int = 1) -> chargefw.Molecule:
    coordinates = [
        [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
        [[0.1, 0.0, 0.0], [1.06, 0.0, 0.0], [-0.14, 0.93, 0.0]],
    ][:conformers]
    return chargefw.Molecule(
        [8, 1, 1],
        bonds=[[0, 1, 1], [0, 2, 1]],
        coordinates=coordinates,
        name="water",
        source_name="water-input",
        record_id="water-1",
    )


class GeneratedOutputTests(unittest.TestCase):
    def test_generated_formats_use_native_writers(self) -> None:
        molecule = water()
        result = chargefw.calculate(molecule, method="formal")

        sdf = chargefw.io.dumps(result, format="sdf")
        self.assertIn("V3000", sdf)
        self.assertIn("CHARGEFW_CHARGES_1", sdf)
        self.assertIn("method=formal", sdf)

        mol2 = chargefw.io.dumps(result, format="mol2")
        self.assertIn("@<TRIPOS>MOLECULE\nwater", mol2)
        self.assertIn("USER_CHARGES", mol2)

        mmcif = chargefw.io.dumps(result, format="mmcif")
        self.assertIn("data_water-1", mmcif)
        self.assertIn("_sb_ncbr_partial_atomic_charges.", mmcif)

        encoded = json.loads(chargefw.io.dumps(result, format="result-json"))
        self.assertEqual(encoded["status"], "success")
        self.assertEqual(encoded["results"][0]["input"]["source"], "water-input")
        self.assertEqual(encoded["results"][0]["input"]["record_id"], "water-1")
        self.assertEqual(encoded["calculation_provenance"]["effective"]["method"]["id"], "formal")
        self.assertIs(result.molecules[0], molecule)

    def test_write_uses_explicit_format(self) -> None:
        result = chargefw.calculate(water(), method="formal")
        with TemporaryDirectory() as directory:
            path = Path(directory) / "charges.data"
            chargefw.io.write(path, result, format="mol2")
            self.assertIn("@<TRIPOS>MOLECULE", path.read_text(encoding="utf-8"))

    def test_single_geometry_output_requires_unambiguous_conformer(self) -> None:
        result = chargefw.calculate(water(conformers=2), method="formal")
        with self.assertRaisesRegex(ValueError, "conformer is required"):
            chargefw.io.dumps(result, format="sdf")
        selected = chargefw.io.dumps(result, format="sdf", conformer=1)
        self.assertIn("M  V30 1 O 0.1 0 0", selected)
        self.assertIn("_atom_site.pdbx_PDB_model_num", chargefw.io.dumps(result, format="mmcif"))

    def test_mmcif_applies_geometry_independent_charges_to_selected_conformers(self) -> None:
        contents = json.dumps(
            {
                "schema_version": "1.0",
                "molecules": [
                    {
                        "atoms": [{"atomic_number": 1, "formal_charge": 0}],
                        "conformers": [
                            {"coordinates": [[0, 0, 0]]},
                            {"coordinates": [[1, 0, 0]]},
                        ],
                    }
                ],
            }
        )
        for selection, expected_count in (("all", 2), ("first", 1)):
            with self.subTest(selection=selection):
                molecules = chargefw.io.parse(
                    contents, format="molecule-json", conformers=cast(Any, selection)
                )
                result = chargefw.calculate(molecules, method="formal")
                block = gemmi.cif.read_string(
                    chargefw.io.dumps(result, format="mmcif")
                ).sole_block()
                atom_sites = block.find("_atom_site.", ["id"])
                charges = block.find("_sb_ncbr_partial_atomic_charges.", ["atom_id"])
                self.assertEqual(len(atom_sites), expected_count)
                self.assertEqual(
                    [gemmi.cif.as_string(row[0]) for row in charges],
                    [gemmi.cif.as_string(row[0]) for row in atom_sites],
                )

    def test_result_json_supports_failed_calculations(self) -> None:
        molecule = chargefw.Molecule([8, 1, 1], bonds=[[0, 1, 1], [0, 2, 1]])
        try:
            chargefw.calculate(molecule, method="qeq")
        except chargefw.NoExecutablePlanError as error:
            failed_result = error.result
            encoded = json.loads(chargefw.io.dumps(failed_result, format="result-json"))
        else:
            self.fail("calculation unexpectedly succeeded")

        self.assertEqual(encoded["status"], "no_executable_plan")
        self.assertNotIn("assignments", encoded["results"][0])
        self.assertTrue(encoded["results"][0]["diagnostics"])
        with self.assertRaisesRegex(ValueError, "successful calculation"):
            chargefw.io.dumps(failed_result, format="mmcif")

    def test_result_json_preserves_import_diagnostics(self) -> None:
        path = Path(__file__).parents[1] / "fixtures" / "synthetic" / "mol2" / "aromatic.mol2"
        result = chargefw.calculate(chargefw.io.read(path, format="mol2"), method="formal")

        encoded = json.loads(chargefw.io.dumps(result, format="result-json"))

        diagnostics = encoded["results"][0]["diagnostics"]
        self.assertEqual([value["code"] for value in diagnostics], ["partial_charges_ignored"])

    def test_molecular_output_requires_finite_coordinates(self) -> None:
        missing = chargefw.calculate(chargefw.Molecule([1]), method="formal")
        with self.assertRaisesRegex(ValueError, "conformer|coordinates"):
            chargefw.io.dumps(missing, format="sdf")
        with self.assertRaisesRegex(ValueError, "coordinates"):
            chargefw.io.dumps(missing, format="mmcif")

        nonfinite = chargefw.calculate(
            chargefw.Molecule([1], coordinates=[[float("nan"), 0.0, 0.0]]),
            method="formal",
        )
        with self.assertRaisesRegex(ValueError, "must be finite"):
            chargefw.io.dumps(nonfinite, format="mol2")

    def test_output_arguments_are_explicit(self) -> None:
        result = chargefw.calculate(water(), method="formal")
        with self.assertRaises(TypeError):
            chargefw.io.dumps(result)  # type: ignore[call-arg]
        with self.assertRaisesRegex(ValueError, "unsupported calculation output format"):
            chargefw.io.dumps(result, format=cast(Any, "pdb"))
        with self.assertRaisesRegex(ValueError, "only supported for SDF"):
            chargefw.io.dumps(result, format="mol2", sdf_version="v2000")
        with self.assertRaisesRegex(ValueError, "only supported for SDF and MOL2"):
            chargefw.io.dumps(result, format="mmcif", conformer=0)

    def test_manual_molecules_do_not_claim_import_provenance(self) -> None:
        result = chargefw.calculate(water(), method="formal")

        encoded = json.loads(chargefw.io.dumps(result, format="result-json"))

        requested = encoded["calculation_provenance"]["requested"]
        self.assertNotIn("input", requested)
        self.assertNotIn("structural_input", requested)


if __name__ == "__main__":
    unittest.main()
