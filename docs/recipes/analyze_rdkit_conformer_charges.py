"""Analyze calculated charge variation across an RDKit conformer set.

Run with: python docs/recipes/analyze_rdkit_conformer_charges.py "CCO" --conformers 20
"""

from __future__ import annotations

import argparse
from typing import Any

import chargefw
import numpy as np
from chargefw.io import rdkit as chargefw_rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


def prepare_demo_ensemble(smiles: str, conformer_count: int) -> Any:
    """Generate an unoptimized, deterministic ETKDG ensemble for demonstration."""

    if conformer_count < 2:
        raise ValueError("conformer_count must be at least 2")
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")

    Chem.Kekulize(molecule, clearAromaticFlags=True)
    molecule = Chem.AddHs(molecule)
    embedding = AllChem.ETKDGv3()
    embedding.randomSeed = 0xC0FFEE
    conformer_ids = AllChem.EmbedMultipleConfs(molecule, numConfs=conformer_count, params=embedding)
    if len(conformer_ids) < 2:
        raise RuntimeError(f"RDKit generated only {len(conformer_ids)} conformers")
    return molecule


def atom_charge_statistics(
    molecule: Any, *, source_name: str = ""
) -> tuple[list[tuple[str, float, float, float, float]], chargefw.CalculationResult]:
    """Calculate per-atom mean, standard deviation, minimum, and maximum charge."""

    native_molecule = chargefw_rdkit.from_mol(molecule, source_name=source_name)
    result = chargefw.calculate(
        native_molecule,
        method="qeq",
        parameter_set="QEq_original",
        options={"overlap_term": "Louwen-Vogt"},
        execution="full",
    )
    charges = np.stack([assignment.values for assignment in result.assignments])
    statistics = []
    for atom in molecule.GetAtoms():
        values = charges[:, atom.GetIdx()]
        statistics.append(
            (
                f"{atom.GetSymbol()}{atom.GetIdx() + 1}",
                float(np.mean(values)),
                float(np.std(values)),
                float(np.min(values)),
                float(np.max(values)),
            )
        )
    return statistics, result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smiles", nargs="?", default="CCO")
    parser.add_argument("--conformers", type=int, default=20)
    parser.add_argument("--top", type=int, default=5)
    arguments = parser.parse_args()

    molecule = prepare_demo_ensemble(arguments.smiles, arguments.conformers)
    statistics, result = atom_charge_statistics(molecule, source_name=arguments.smiles)
    parameter_set = result.plan.parameter_set.id if result.plan.parameter_set else "-"
    print(f"Method: {result.plan.method.id}/{parameter_set}")
    print(f"Conformers: {len(result.assignments)} equally weighted, unoptimized ETKDGv3 structures")
    print("Atom          mean         std         min         max")
    for label, mean, standard_deviation, minimum, maximum in sorted(
        statistics, key=lambda item: item[2], reverse=True
    )[: arguments.top]:
        print(f"{label:<6} {mean:12.6f} {standard_deviation:11.6f} {minimum:11.6f} {maximum:11.6f}")


if __name__ == "__main__":
    main()
