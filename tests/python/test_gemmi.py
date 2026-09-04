"""Molecular format input and Gemmi conversion checks."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, cast

import gemmi
import numpy as np
from chargefw import io as chargefw_io
from chargefw.io import gemmi as chargefw_gemmi

PDB_TEXT = """HEADER    TEST PDB
TITLE     TWO MODELS
MODEL        1
ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O  
ATOM      2  H1 AHOH A   1       0.957   0.000   0.000  1.00 20.00           H  
ATOM      3  H1 BHOH A   1       9.000   0.000   0.000  1.00 20.00           H  
ENDMDL
MODEL        2
ATOM      1  O   HOH A   1       0.100   0.000   0.000  1.00 20.00           O  
ATOM      2  H1 AHOH A   1       1.057   0.000   0.000  1.00 20.00           H  
ATOM      3  H1 BHOH A   1       9.100   0.000   0.000  1.00 20.00           H  
ENDMDL
CONECT    1    2
END
"""

MMCIF_TEXT = """data_first
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 ALA A CA 1
HETATM 2 O O . HOH A 2 ? 1.0 0.0 0.0 1.0 20.0 0 2 HOH A O 1
ATOM 3 C CA . ALA A 1 ? 0.1 0.0 0.0 1.0 20.0 0 1 ALA A CA 2
HETATM 4 O O . HOH A 2 ? 1.1 0.0 0.0 1.0 20.0 0 2 HOH A O 2
#
data_second
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 O O . HOH B 1 ? 3.0 0.0 0.0 1.0 20.0 0 1 HOH B O 1
#
"""

BOND_STRATEGY_PDB = """SSBOND   1 CYS A   3    CYS A   4                          
LINK         C   ALA A   1                 C1  LIG A   5
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N  
ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C  
ATOM      3  C   ALA A   1       2.450   1.000   0.000  1.00 20.00           C  
ATOM      4  O   ALA A   1       3.450   1.000   0.000  1.00 20.00           O  
ATOM      5  CB  ALA A   1       1.450  -1.000   0.000  1.00 20.00           C  
ATOM      6  N   GLY A   2       2.200   2.200   0.000  1.00 20.00           N  
ATOM      7  CA  GLY A   2       3.200   3.200   0.000  1.00 20.00           C  
ATOM      8  C   GLY A   2       4.200   3.200   0.000  1.00 20.00           C  
ATOM      9  O   GLY A   2       5.200   3.200   0.000  1.00 20.00           O  
ATOM     10  SG  CYS A   3       6.200   3.200   0.000  1.00 20.00           S  
ATOM     11  SG  CYS A   4       7.200   3.200   0.000  1.00 20.00           S  
HETATM   12  C1  LIG A   5       8.200   3.200   0.000  1.00 20.00           C  
CONECT   11   12
END
"""

MOL_TEXT = """charged-v2000
  ChargeFW

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  CHG  2   1   1   2  -1
M  END
"""

MOL2_TEXT = """@<TRIPOS>MOLECULE
water
3 2 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
1 O 0.0000 0.0000 0.0000 O.2
2 H1 0.9570 0.0000 0.0000 H
3 H2 -0.2400 0.9270 0.0000 H
@<TRIPOS>BOND
1 1 2 1
2 1 3 1
"""

MOLECULE_JSON_TEXT = """{
  "schema_version": "1.0",
  "molecules": [{
    "id": "water-1",
    "name": "water",
    "atoms": [
      {"atomic_number": 8, "formal_charge": 0},
      {"atomic_number": 1, "formal_charge": 0}
    ],
    "bonds": [{"atoms": [0, 1], "order": 1}],
    "conformers": [
      {"id": "first", "coordinates": [[0, 0, 0], [1, 0, 0]]},
      {"id": "second", "coordinates": [[0.1, 0, 0], [1.1, 0, 0]]}
    ]
  }]
}
"""


class NativeInputTests(unittest.TestCase):
    def test_parse_native_molecular_formats(self) -> None:
        molecule = chargefw_io.parse_mol(MOL_TEXT, source_name="charged.mol")
        self.assertEqual(molecule.atomic_numbers.tolist(), [7, 8])
        self.assertEqual(molecule.formal_charges.tolist(), [1, -1])
        self.assertEqual(molecule.source_name, "charged.mol")

        sdf = chargefw_io.parse_sdf(f"{MOL_TEXT}$$$$\n", source_name="charged.sdf")
        self.assertEqual(len(sdf), 1)
        self.assertEqual(sdf[0].record_index, 0)

        mol2 = chargefw_io.parse_mol2(MOL2_TEXT, source_name="water.mol2")
        self.assertEqual(len(mol2), 1)
        self.assertEqual(mol2[0].atom_names, ("O", "H1", "H2"))
        self.assertEqual(mol2[0].bond_count, 2)

    def test_parse_molecule_json_conformer_selection(self) -> None:
        all_conformers = chargefw_io.parse_molecule_json(MOLECULE_JSON_TEXT)
        first_conformer = chargefw_io.parse_molecule_json(
            MOLECULE_JSON_TEXT, conformers="first"
        )

        self.assertEqual(all_conformers[0].conformer_names, ("first", "second"))
        self.assertEqual(first_conformer[0].conformer_names, ("first",))

    def test_path_readers_set_source_name(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "charged.mol"
            path.write_text(MOL_TEXT, encoding="utf-8")
            molecule = chargefw_io.read_mol(path)

        self.assertEqual(molecule.source_name, str(path))


class GemmiAdapterTests(unittest.TestCase):
    def test_pdb_text_file_and_structure_preserve_mapping(self) -> None:
        text_molecule = chargefw_io.parse_pdb(
            PDB_TEXT,
            source_name="water.pdb",
            bonds="explicit",
        )
        self.assertEqual(text_molecule.atom_names, ("O", "H1"))
        self.assertEqual(text_molecule.atom_ids, (0, 1))
        self.assertEqual(text_molecule.conformer_names, ("1", "2"))
        self.assertEqual(text_molecule.source_name, "water.pdb")
        self.assertEqual(text_molecule.bonds.tolist(), [[0, 1, 1]])
        np.testing.assert_allclose(text_molecule.coordinates[:, 0, 0], [0.0, 0.1])

        with TemporaryDirectory() as directory:
            path = Path(directory) / "water.pdb"
            path.write_text(PDB_TEXT, encoding="utf-8")
            file_molecule = chargefw_io.read_pdb(path, bonds="explicit")
        np.testing.assert_array_equal(file_molecule.atomic_numbers, text_molecule.atomic_numbers)
        np.testing.assert_array_equal(file_molecule.coordinates, text_molecule.coordinates)

        structure = gemmi.read_pdb_string(PDB_TEXT)
        structure_molecule = chargefw_gemmi.from_structure(
            structure, source_name="structure"
        )
        np.testing.assert_array_equal(
            structure_molecule.atomic_numbers, text_molecule.atomic_numbers
        )
        np.testing.assert_array_equal(structure_molecule.coordinates, text_molecule.coordinates)

    def test_mmcif_text_file_and_document_preserve_records(self) -> None:
        collection = chargefw_io.parse_mmcif(MMCIF_TEXT, source_name="models.cif")
        self.assertEqual(len(collection), 2)
        self.assertEqual([value.record_index for value in collection], [0, 1])
        self.assertEqual([value.record_id for value in collection], ["first", "second"])
        self.assertEqual(collection[0].conformer_names, ("1", "2"))

        with TemporaryDirectory() as directory:
            path = Path(directory) / "models.cif"
            path.write_text(MMCIF_TEXT, encoding="utf-8")
            from_file = chargefw_io.read_mmcif(path)
        self.assertEqual([value.atom_count for value in from_file], [2, 1])

        document = gemmi.cif.read_string(MMCIF_TEXT)
        from_document = chargefw_gemmi.from_document(document, source_name="document")
        self.assertEqual([value.record_id for value in from_document], ["first", "second"])
        np.testing.assert_array_equal(from_document[0].coordinates, collection[0].coordinates)

    def test_selection_conformers_and_types_are_explicit(self) -> None:
        polymers = chargefw_io.parse_mmcif(
            MMCIF_TEXT,
            selection="polymers",
            conformers="first",
        )
        self.assertEqual(polymers[0].atom_names, ("CA",))
        self.assertEqual(polymers[0].conformer_count, 1)

        with self.assertRaises(TypeError):
            chargefw_gemmi.from_structure(cast(Any, object()))
        with self.assertRaises(ValueError):
            chargefw_io.parse_pdb(PDB_TEXT, selection=cast(Any, "invalid"))
        with self.assertRaises(ValueError):
            chargefw_io.parse_pdb(
                PDB_TEXT, selection=cast(Any, "polymers_and_ligands")
            )
        with self.assertRaises(TypeError):
            chargefw_io.parse_pdb(PDB_TEXT, selection=cast(Any, 1))

    def test_bond_strategies_match_native_adapter(self) -> None:
        expected_counts = {
            "none": 0,
            "templates": 8,
            "explicit": 2,
            "hybrid": 10,
        }
        for strategy, expected_count in expected_counts.items():
            with self.subTest(strategy=strategy):
                molecule = chargefw_io.parse_pdb(
                    BOND_STRATEGY_PDB, bonds=cast(Any, strategy)
                )
                self.assertEqual(molecule.bond_count, expected_count)

    def test_source_atom_ids_distinguish_repeated_atom_names(self) -> None:
        molecule = chargefw_io.parse_pdb(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
            "ATOM      2  CA  GLY A   2       1.000   0.000   0.000  1.00 20.00           C  \n"
            "END\n"
        )
        self.assertEqual(molecule.atom_names, ("CA", "CA"))
        self.assertEqual(molecule.atom_ids, (0, 1))


if __name__ == "__main__":
    unittest.main()
