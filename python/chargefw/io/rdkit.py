"""Optional conversion and charge attachment for RDKit molecules."""

from __future__ import annotations

from importlib import import_module
from operator import index as as_index
from typing import TYPE_CHECKING, Any

from ..core import Molecule

if TYPE_CHECKING:
    from ..calculation import CalculationResult
    from ..charges import ChargeAssignment


def _require_rdkit() -> Any:
    try:
        return import_module("rdkit.Chem")
    except ModuleNotFoundError as error:
        if error.name not in ("rdkit", "rdkit.Chem"):
            raise
        raise ImportError(
            "chargefw.io.rdkit requires an independently installed RDKit package"
        ) from error


def from_mol(molecule: Any, *, source_name: str = "") -> Molecule:
    """Copy an existing RDKit molecule without preparation or sanitization."""

    chemistry = _require_rdkit()
    if not isinstance(molecule, chemistry.Mol):
        raise TypeError("molecule must be an rdkit.Chem.Mol")
    if not isinstance(source_name, str):
        raise TypeError("source_name must be a string")

    atoms = tuple(molecule.GetAtoms())
    bonds: list[tuple[int, int, int]] = []
    for bond in molecule.GetBonds():
        order_value = bond.GetBondTypeAsDouble()
        if order_value not in (1.0, 2.0, 3.0):
            raise ValueError(
                "RDKit molecule contains a bond that is not explicitly single, double, or triple"
            )
        bonds.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), int(order_value)))

    coordinates = [
        [
            (
                conformer.GetAtomPosition(index).x,
                conformer.GetAtomPosition(index).y,
                conformer.GetAtomPosition(index).z,
            )
            for index in range(len(atoms))
        ]
        for conformer in molecule.GetConformers()
    ]
    atom_names = tuple(
        atom.GetProp("_TriposAtomName")
        if atom.HasProp("_TriposAtomName")
        else f"{atom.GetSymbol()}{atom.GetIdx() + 1}"
        for atom in atoms
    )
    name = molecule.GetProp("_Name") if molecule.HasProp("_Name") else ""
    return Molecule(
        [atom.GetAtomicNum() for atom in atoms],
        formal_charges=[atom.GetFormalCharge() for atom in atoms],
        bonds=bonds,
        coordinates=coordinates,
        name=name,
        atom_names=atom_names,
        conformer_names=[str(value.GetId()) for value in molecule.GetConformers()],
        source_name=source_name,
        atom_ids=range(len(atoms)),
    )


def _assignment(
    result: CalculationResult, molecule_index: int, conformer: int | None
) -> ChargeAssignment:
    assignments = result.assignments_by_molecule[molecule_index]
    if conformer is not None:
        matches = tuple(value for value in assignments if value.conformer_index == conformer)
        if len(matches) == 1:
            return matches[0]
        if len(assignments) == 1 and assignments[0].conformer_index is None:
            return assignments[0]
        raise ValueError("selected conformer has no charge assignment")
    if len(assignments) != 1:
        raise ValueError("conformer is required when the result has multiple assignments")
    return assignments[0]


def attach_charges(
    molecule: Any,
    result: CalculationResult,
    *,
    molecule_index: int = 0,
    conformer: int | None = None,
    property_name: str = "ChargeFWPartialCharge",
    overwrite: bool = False,
) -> None:
    """Attach one calculated charge assignment to RDKit atom properties in place."""

    from ..calculation import CalculationResult

    chemistry = _require_rdkit()
    if not isinstance(molecule, chemistry.Mol):
        raise TypeError("molecule must be an rdkit.Chem.Mol")
    if not isinstance(result, CalculationResult):
        raise TypeError("result must be a CalculationResult")
    if isinstance(molecule_index, bool):
        raise TypeError("molecule_index must be an integer")
    try:
        molecule_index = as_index(molecule_index)
    except TypeError as error:
        raise TypeError("molecule_index must be an integer") from error
    if molecule_index < 0 or molecule_index >= len(result.molecules):
        raise IndexError("molecule_index is outside the calculation input collection")
    if isinstance(conformer, bool):
        raise TypeError("conformer must be an integer or None")
    if conformer is not None:
        try:
            conformer = as_index(conformer)
        except TypeError as error:
            raise TypeError("conformer must be an integer or None") from error
        if conformer < 0:
            raise ValueError("conformer must be non-negative")
    if not isinstance(property_name, str):
        raise TypeError("property_name must be a string")
    if not property_name:
        raise ValueError("property_name must not be empty")
    if not isinstance(overwrite, bool):
        raise TypeError("overwrite must be a bool")

    assignment = _assignment(result, molecule_index, conformer)
    source = result.molecules[molecule_index]
    if molecule.GetNumAtoms() != source.atom_count:
        raise ValueError("RDKit target atom count does not match the calculation input")
    property_list_name = f"atom.dprop.{property_name}"
    if not overwrite and molecule.HasProp(property_list_name):
        raise ValueError(f"RDKit molecule already has property {property_list_name!r}")
    for source_index, atom_id in enumerate(source.atom_ids):
        if not isinstance(atom_id, int) or isinstance(atom_id, bool):
            raise ValueError("RDKit attachment requires integer atom IDs")
        if atom_id < 0 or atom_id >= molecule.GetNumAtoms():
            raise ValueError("RDKit atom ID is outside the target molecule")
        atom = molecule.GetAtomWithIdx(atom_id)
        if atom.GetAtomicNum() != int(source.atomic_numbers[source_index]):
            raise ValueError("RDKit target atom mapping does not match atomic numbers")
        if atom.GetFormalCharge() != int(source.formal_charges[source_index]):
            raise ValueError("RDKit target atom mapping does not match formal charges")
        if not overwrite and atom.HasProp(property_name):
            raise ValueError(f"RDKit atom already has property {property_name!r}")

    for atom_id, charge in zip(source.atom_ids, assignment.values, strict=True):
        molecule.GetAtomWithIdx(atom_id).SetDoubleProp(property_name, float(charge))
    create_property_list = getattr(chemistry, "CreateAtomDoublePropertyList", None)
    if create_property_list is not None:
        create_property_list(molecule, property_name)


__all__ = ["from_mol", "attach_charges"]
