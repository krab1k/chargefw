"""Find atoms whose calculated charge changes most across RDKit conformers.

Run with: python docs/recipes/conformational_charge_hotspots.py "CCO" --conformers 20
"""

from __future__ import annotations

import argparse

import chargefw
import numpy as np
from chargefw.io import rdkit as chargefw_rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


def atom_charge_variances(smiles: str, conformer_count: int = 20) -> list[tuple[str, float]]:
    """Generate conformers and return each atom label with its charge variance."""

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    if conformer_count < 2:
        raise ValueError("conformer_count must be at least 2")

    Chem.Kekulize(molecule, clearAromaticFlags=True)
    molecule = Chem.AddHs(molecule)
    embedding = AllChem.ETKDGv3()
    embedding.randomSeed = 0xC0FFEE
    conformer_ids = AllChem.EmbedMultipleConfs(
        molecule, numConfs=conformer_count, params=embedding
    )
    if len(conformer_ids) != conformer_count:
        raise RuntimeError(f"RDKit generated only {len(conformer_ids)} conformers")

    native_molecule = chargefw_rdkit.from_mol(molecule, source_name=smiles)
    result = chargefw.calculate(native_molecule, method="qeq", execution="full")
    charges = np.stack([assignment.values for assignment in result.assignments])
    variances = np.var(charges, axis=0)
    return [
        (f"{atom.GetSymbol()}{atom.GetIdx() + 1}", float(variances[atom.GetIdx()]))
        for atom in molecule.GetAtoms()
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smiles", nargs="?", default="CCO")
    parser.add_argument("--conformers", type=int, default=20)
    arguments = parser.parse_args()

    variances = atom_charge_variances(arguments.smiles, arguments.conformers)
    lowest = min(variances, key=lambda item: item[1])
    highest = max(variances, key=lambda item: item[1])
    print(f"Lowest charge variance:  {lowest[0]} ({lowest[1]:.6g})")
    print(f"Highest charge variance: {highest[0]} ({highest[1]:.6g})")


if __name__ == "__main__":
    main()
