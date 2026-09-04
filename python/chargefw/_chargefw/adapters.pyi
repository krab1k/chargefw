from typing import Literal, TypedDict

class MoleculePayload(TypedDict):
    atomic_numbers: list[int]
    formal_charges: list[int]
    bonds: list[tuple[int, int, int]]
    coordinates: list[list[tuple[float, float, float]]]
    name: str
    atom_names: list[str]
    conformer_names: list[str]
    source: str
    record_index: int
    record_id: str

def _parse_mol(contents: str, source: str) -> list[MoleculePayload]: ...
def _parse_sdf(contents: str, source: str) -> list[MoleculePayload]: ...
def _parse_mol2(contents: str, source: str) -> list[MoleculePayload]: ...
def _parse_molecule_json(
    contents: str,
    source: str,
    conformers: Literal["first", "all"],
) -> list[MoleculePayload]: ...
def _parse_pdb(
    contents: str,
    source: str,
    selection: Literal["all", "polymers-and-ligands", "polymers"],
    bonds: Literal["none", "explicit", "templates", "hybrid"],
    conformers: Literal["first", "all"],
) -> list[MoleculePayload]: ...
def _parse_mmcif(
    contents: str,
    source: str,
    selection: Literal["all", "polymers-and-ligands", "polymers"],
    bonds: Literal["none", "explicit", "templates", "hybrid"],
    conformers: Literal["first", "all"],
) -> list[MoleculePayload]: ...
