from enum import Enum
from typing import TypedDict

class RecordSelection(Enum):
    ALL = ...
    POLYMERS_AND_LIGANDS = ...
    POLYMERS = ...


class BondStrategy(Enum):
    NONE = ...
    EXPLICIT_BONDS = ...
    TEMPLATES = ...
    HYBRID = ...


class ConformerSelection(Enum):
    FIRST = ...
    ALL = ...


class DiagnosticPayload(TypedDict):
    code: str
    message: str
    line: int | None


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
    diagnostics: list[DiagnosticPayload]


def _read_pdb(
    contents: str,
    source: str,
    selection: RecordSelection,
    bond_strategy: BondStrategy,
    conformers: ConformerSelection,
) -> list[MoleculePayload]: ...


def _read_mmcif(
    contents: str,
    source: str,
    selection: RecordSelection,
    bond_strategy: BondStrategy,
    conformers: ConformerSelection,
) -> list[MoleculePayload]: ...
