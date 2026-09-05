from collections.abc import Sequence
from typing import Literal, TypedDict

from .calculation import _NativeExecutionResult
from .core import _NativeMolecule

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
    diagnostics: list[tuple[str, str, int | None]]

def _parse(
    contents: str,
    source: str,
    format: Literal["mol", "sdf", "mol2", "molecule-json", "pdb", "mmcif"],
    selection: Literal["all", "polymers-and-ligands", "polymers"],
    bonds: Literal["none", "explicit", "templates", "hybrid"],
    conformers: Literal["first", "all"],
) -> list[MoleculePayload]: ...
def _dumps(
    result: _NativeExecutionResult,
    molecules: Sequence[_NativeMolecule],
    identities: Sequence[tuple[str, int, str]],
    diagnostics: Sequence[Sequence[tuple[str, str, int | None]]],
    requested: dict[str, object],
    format: Literal["sdf", "mol2", "mmcif", "result-json"],
    conformer: int | None,
    sdf_version: Literal["v2000", "v3000"],
) -> str: ...
def _attach_mmcif(
    contents: str,
    result: _NativeExecutionResult,
    molecules: Sequence[_NativeMolecule],
    selection: Literal["all", "polymers-and-ligands", "polymers"],
    overwrite: bool,
) -> str: ...
