"""Inspect the normalized molecules imported by ChargeFW.

Run with: python docs/recipes/inspect_molecules.py input.sdf --format sdf
"""

from __future__ import annotations

import argparse
from pathlib import Path

import chargefw


def inspect_molecules(
    input_path: str | Path,
    *,
    input_format: str,
    selection: str = "all",
    bonds: str = "none",
    conformers: str = "all",
) -> None:
    """Print source identities, dimensions, and formal charges for imported records."""

    molecules = chargefw.io.read(
        input_path,
        format=input_format,
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )
    print(f"Collection: {molecules.name or '<unnamed>'} ({len(molecules)} molecule(s))")

    for index, molecule in enumerate(molecules):
        print(
            f"[{index}] record={molecule.record_index} id={molecule.record_id!r} "
            f"name={molecule.name!r}"
        )
        print(
            f"    atoms={molecule.atom_count} bonds={molecule.bond_count} "
            f"conformers={molecule.conformer_count} "
            f"formal_charge={int(molecule.formal_charges.sum())}"
        )
        print(f"    conformer_names={molecule.conformer_names!r}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--format", choices=chargefw.io.INPUT_FORMATS, required=True)
    parser.add_argument(
        "--selection", choices=("all", "polymers-and-ligands", "polymers"), default="all"
    )
    parser.add_argument(
        "--bonds", choices=("none", "explicit", "templates", "hybrid"), default="none"
    )
    parser.add_argument("--conformers", choices=("first", "all"), default="all")
    arguments = parser.parse_args()
    inspect_molecules(
        arguments.input,
        input_format=arguments.format,
        selection=arguments.selection,
        bonds=arguments.bonds,
        conformers=arguments.conformers,
    )


if __name__ == "__main__":
    main()
