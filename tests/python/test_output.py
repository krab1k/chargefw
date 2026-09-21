"""Generated molecular and calculation-result output checks."""

import json
import unittest
from gc import collect
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, cast

import chargefw
import gemmi
import numpy as np


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


class OutputTests(unittest.TestCase):
    def test_formats_use_native_writers(self) -> None:
        molecule = water()
        result = chargefw.calculate(molecule, method="formal")

        mol2 = chargefw.io.dumps(result, format="mol2")
        self.assertIn("@<TRIPOS>MOLECULE\nwater-1", mol2)
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
        molecule = chargefw.io.parse(
            """{
  "schema_version": "1.0",
  "molecules": [{
    "atoms": [
      {"atomic_number": 8, "formal_charge": 0},
      {"atomic_number": 1, "formal_charge": 0},
      {"atomic_number": 1, "formal_charge": 0}
    ],
    "bonds": [
      {"atoms": [0, 1], "order": 1},
      {"atoms": [0, 2], "order": 1}
    ]
  }]
}""",
            format="molecule-json",
            source_name="topology.json",
        )[0]
        try:
            chargefw.calculate(molecule, method="qeq")
        except chargefw.NoExecutablePlanError as error:
            failed_result = error.result
            encoded = json.loads(chargefw.io.dumps(failed_result, format="result-json"))
        else:
            self.fail("calculation unexpectedly succeeded")

        self.assertEqual(encoded["status"], "no_executable_plan")
        self.assertNotIn("assignments", encoded["results"][0])
        imported = encoded["results"][0]["input"]["import"]
        self.assertEqual(imported["format"], "molecule-json")
        self.assertEqual(len(imported["atom_mapping"]), 3)
        self.assertEqual(imported["conformer_mapping"], [])
        self.assertTrue(encoded["results"][0]["diagnostics"])
        for format_name in ("mol2", "mmcif"):
            with (
                self.subTest(format=format_name),
                self.assertRaisesRegex(ValueError, "successful calculation"),
            ):
                chargefw.io.dumps(failed_result, format=cast(Any, format_name))

    def test_result_json_preserves_import_diagnostics(self) -> None:
        path = Path(__file__).parents[1] / "fixtures" / "synthetic" / "mol2" / "aromatic.mol2"
        result = chargefw.calculate(chargefw.io.read(path, format="mol2"), method="formal")

        encoded = json.loads(chargefw.io.dumps(result, format="result-json"))

        diagnostics = encoded["results"][0]["diagnostics"]
        self.assertEqual([value["code"] for value in diagnostics], ["partial_charges_ignored"])

    def test_import_metadata_follows_reordered_and_recombined_molecules(self) -> None:
        path = Path(__file__).parents[1] / "fixtures" / "synthetic" / "mol2" / "aromatic.mol2"
        imported = chargefw.io.read(path, format="mol2")
        imported_molecule = imported[0]
        del imported
        collect()

        manual = chargefw.Molecule([1], source_name="manual")
        result = chargefw.calculate([manual, imported_molecule], method="formal")
        encoded = json.loads(chargefw.io.dumps(result, format="result-json"))

        self.assertEqual(encoded["results"][0]["input"]["source"], "manual")
        self.assertNotIn("import", encoded["results"][0]["input"])
        self.assertEqual(encoded["results"][0]["diagnostics"], [])
        imported_input = encoded["results"][1]["input"]["import"]
        self.assertEqual(imported_input["format"], "mol2")
        self.assertEqual(
            [value["source_id"] for value in imported_input["atom_mapping"]],
            list(imported_molecule.atom_ids),
        )
        self.assertEqual(
            [value["code"] for value in encoded["results"][1]["diagnostics"]],
            ["partial_charges_ignored"],
        )
        requested = encoded["calculation_provenance"]["requested"]
        self.assertNotIn("input", requested)
        self.assertNotIn("structural_input", requested)

    def test_coordinate_free_result_json_retains_mapping(self) -> None:
        molecules = chargefw.io.parse(
            """{
  "schema_version": "1.0",
  "molecules": [{"atoms": [{"atomic_number": 1, "formal_charge": 0}]}]
}""",
            format="molecule-json",
            source_name="hydrogen.json",
        )
        result = chargefw.calculate(molecules, method="formal")
        del molecules
        collect()

        encoded = json.loads(chargefw.io.dumps(result, format="result-json"))
        record = encoded["results"][0]
        self.assertEqual(record["assignments"][0]["scope"], "molecule")
        self.assertEqual(record["input"]["import"]["atom_mapping"], [{"source_position": 0}])
        self.assertEqual(record["input"]["import"]["conformer_mapping"], [])

    def test_result_json_preserves_integer_caller_ids(self) -> None:
        molecule = chargefw.Molecule(
            [1, 1], record_id=42, atom_ids=cast(Any, ["left", np.int64(9)])
        )
        result = chargefw.calculate(molecule, method="formal")

        encoded = json.loads(chargefw.io.dumps(result, format="result-json"))
        self.assertEqual(encoded["results"][0]["input"]["record_id"], 42)
        self.assertEqual(encoded["results"][0]["input"]["atom_ids"], ["left", 9])
        self.assertEqual(result.assignments[0].atom_ids, ("left", 9))

    def test_failed_result_json_preserves_caller_atom_ids(self) -> None:
        molecule = chargefw.Molecule([1], atom_ids=["hydrogen"])
        try:
            chargefw.calculate(molecule, method="qeq")
        except chargefw.NoExecutablePlanError as error:
            encoded = json.loads(chargefw.io.dumps(error.result, format="result-json"))
        else:
            self.fail("calculation unexpectedly succeeded")

        self.assertEqual(encoded["results"][0]["input"]["atom_ids"], ["hydrogen"])
        self.assertNotIn("assignments", encoded["results"][0])

    def test_result_json_omits_default_and_imported_atom_ids(self) -> None:
        manual = chargefw.calculate(chargefw.Molecule([1]), method="formal")
        manual_input = json.loads(chargefw.io.dumps(manual, format="result-json"))["results"][0][
            "input"
        ]
        self.assertNotIn("atom_ids", manual_input)

        path = Path(__file__).parents[1] / "fixtures" / "synthetic" / "mol2" / "aromatic.mol2"
        imported = chargefw.calculate(chargefw.io.read(path, format="mol2"), method="formal")
        imported_input = json.loads(chargefw.io.dumps(imported, format="result-json"))["results"][
            0
        ]["input"]
        self.assertNotIn("atom_ids", imported_input)
        self.assertIn("atom_mapping", imported_input["import"])

    def test_molecular_output_requires_finite_coordinates(self) -> None:
        missing = chargefw.calculate(chargefw.Molecule([1]), method="formal")
        result_json = chargefw.io.dumps(missing, format="result-json")
        for format_name in ("mol2", "mmcif"):
            with (
                self.subTest(format=format_name),
                self.assertRaisesRegex(ValueError, "coordinates"),
            ):
                chargefw.io.dumps(missing, format=cast(Any, format_name))
        self.assertEqual(chargefw.io.dumps(missing, format="result-json"), result_json)

        nonfinite = chargefw.calculate(
            chargefw.Molecule([1], coordinates=[[float("nan"), 0.0, 0.0]]),
            method="formal",
        )
        for format_name in ("mol2", "mmcif"):
            with (
                self.subTest(format=format_name),
                self.assertRaisesRegex(ValueError, "must be finite"),
            ):
                chargefw.io.dumps(nonfinite, format=cast(Any, format_name))

    def test_output_format_is_explicit(self) -> None:
        result = chargefw.calculate(water(), method="formal")
        with self.assertRaises(TypeError):
            chargefw.io.dumps(result)  # type: ignore[call-arg]

    def test_manual_molecules_do_not_claim_import_provenance(self) -> None:
        result = chargefw.calculate(water(), method="formal")

        encoded = json.loads(chargefw.io.dumps(result, format="result-json"))

        requested = encoded["calculation_provenance"]["requested"]
        self.assertNotIn("input", requested)
        self.assertNotIn("structural_input", requested)

    def test_integer_record_ids_are_preserved_by_outputs(self) -> None:
        molecule = chargefw.Molecule([1], coordinates=[[0, 0, 0]], name="hydrogen", record_id=7)
        result = chargefw.calculate(molecule, method="formal")

        self.assertIn("\n7\n", chargefw.io.dumps(result, format="mol2"))
        self.assertIn("data_7", chargefw.io.dumps(result, format="mmcif"))
        encoded = json.loads(chargefw.io.dumps(result, format="result-json"))
        self.assertEqual(encoded["results"][0]["input"]["record_id"], 7)


if __name__ == "__main__":
    unittest.main()
