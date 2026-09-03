"""Gemmi serialization adapter parity and mapping checks."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, cast

import gemmi
import numpy as np
from chargefw.adapters import gemmi as chargefw_gemmi

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


class GemmiAdapterTests(unittest.TestCase):
    def test_pdb_text_file_and_structure_preserve_mapping(self) -> None:
        text_molecule = chargefw_gemmi.read_pdb_string(
            PDB_TEXT,
            source="water.pdb",
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
            file_molecule = chargefw_gemmi.read_pdb(path, bonds="explicit")
        np.testing.assert_array_equal(file_molecule.atomic_numbers, text_molecule.atomic_numbers)
        np.testing.assert_array_equal(file_molecule.coordinates, text_molecule.coordinates)

        structure = gemmi.read_pdb_string(PDB_TEXT)
        structure_molecule = chargefw_gemmi.from_structure(structure, source="structure")
        np.testing.assert_array_equal(
            structure_molecule.atomic_numbers, text_molecule.atomic_numbers
        )
        np.testing.assert_array_equal(structure_molecule.coordinates, text_molecule.coordinates)

    def test_mmcif_text_file_and_document_preserve_records(self) -> None:
        collection = chargefw_gemmi.read_mmcif_string(MMCIF_TEXT, source="models.cif")
        self.assertEqual(len(collection), 2)
        self.assertEqual([value.record_index for value in collection], [0, 1])
        self.assertEqual([value.record_id for value in collection], ["first", "second"])
        self.assertEqual(collection[0].conformer_names, ("1", "2"))

        with TemporaryDirectory() as directory:
            path = Path(directory) / "models.cif"
            path.write_text(MMCIF_TEXT, encoding="utf-8")
            from_file = chargefw_gemmi.read_mmcif(path)
        self.assertEqual([value.atom_count for value in from_file], [2, 1])

        document = gemmi.cif.read_string(MMCIF_TEXT)
        from_document = chargefw_gemmi.from_document(document, source="document")
        self.assertEqual([value.record_id for value in from_document], ["first", "second"])
        np.testing.assert_array_equal(from_document[0].coordinates, collection[0].coordinates)

    def test_selection_conformers_and_types_are_explicit(self) -> None:
        polymers = chargefw_gemmi.read_mmcif_string(
            MMCIF_TEXT,
            selection="polymers",
            conformers="first",
        )
        self.assertEqual(polymers[0].atom_names, ("CA",))
        self.assertEqual(polymers[0].conformer_count, 1)

        with self.assertRaises(TypeError):
            chargefw_gemmi.from_structure(cast(Any, object()))
        with self.assertRaises(ValueError):
            chargefw_gemmi.read_pdb_string(PDB_TEXT, selection=cast(Any, "invalid"))
        with self.assertRaises(ValueError):
            chargefw_gemmi.read_pdb_string(PDB_TEXT, selection=cast(Any, "polymers_and_ligands"))
        with self.assertRaises(TypeError):
            chargefw_gemmi.read_pdb_string(PDB_TEXT, selection=cast(Any, 1))

    def test_bond_strategies_match_native_adapter(self) -> None:
        expected_counts = {
            "none": 0,
            "templates": 8,
            "explicit": 2,
            "hybrid": 10,
        }
        for strategy, expected_count in expected_counts.items():
            with self.subTest(strategy=strategy):
                molecule = chargefw_gemmi.read_pdb_string(
                    BOND_STRATEGY_PDB, bonds=cast(Any, strategy)
                )
                self.assertEqual(molecule.bond_count, expected_count)

    def test_source_atom_ids_distinguish_repeated_atom_names(self) -> None:
        molecule = chargefw_gemmi.read_pdb_string(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
            "ATOM      2  CA  GLY A   2       1.000   0.000   0.000  1.00 20.00           C  \n"
            "END\n"
        )
        self.assertEqual(molecule.atom_names, ("CA", "CA"))
        self.assertEqual(molecule.atom_ids, (0, 1))


if __name__ == "__main__":
    unittest.main()
