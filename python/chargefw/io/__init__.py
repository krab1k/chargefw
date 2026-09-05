"""Molecular input backed by ChargeFW's native readers."""

from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import Literal, TypeAlias

from .._chargefw import adapters as _native_adapters
from ..core import Molecule, MoleculeCollection

InputFormat: TypeAlias = Literal["mol", "sdf", "mol2", "molecule-json", "pdb", "mmcif"]
RecordSelection: TypeAlias = Literal["all", "polymers-and-ligands", "polymers"]
BondStrategy: TypeAlias = Literal["none", "explicit", "templates", "hybrid"]
ConformerSelection: TypeAlias = Literal["first", "all"]

_INPUT_FORMATS = frozenset(("mol", "sdf", "mol2", "molecule-json", "pdb", "mmcif"))
_STRUCTURAL_FORMATS = frozenset(("pdb", "mmcif"))
_MULTI_CONFORMER_FORMATS = frozenset(("molecule-json", "pdb", "mmcif"))


def _molecule(payload: _native_adapters.MoleculePayload) -> Molecule:
    return Molecule(
        atomic_numbers=payload["atomic_numbers"],
        formal_charges=payload["formal_charges"],
        bonds=payload["bonds"],
        coordinates=payload["coordinates"],
        name=payload["name"],
        atom_names=payload["atom_names"],
        conformer_names=payload["conformer_names"],
        source_name=payload["source"],
        record_index=payload["record_index"],
        record_id=payload["record_id"],
    )


def _collection(
    payloads: list[_native_adapters.MoleculePayload], source_name: str
) -> MoleculeCollection:
    return MoleculeCollection((_molecule(payload) for payload in payloads), name=source_name)


def _validate_options(
    format: InputFormat,
    selection: RecordSelection,
    bonds: BondStrategy,
    conformers: ConformerSelection,
) -> None:
    if not isinstance(format, str):
        raise TypeError("format must be a string")
    if format not in _INPUT_FORMATS:
        raise ValueError(f"unsupported molecular input format: {format}")
    if format not in _STRUCTURAL_FORMATS:
        if selection != "all":
            raise ValueError("selection is only supported for PDB and mmCIF input")
        if bonds != "none":
            raise ValueError("bonds is only supported for PDB and mmCIF input")
    if format not in _MULTI_CONFORMER_FORMATS and conformers != "all":
        raise ValueError("conformers is only supported for molecule JSON, PDB, and mmCIF input")


def parse(
    contents: str,
    *,
    format: InputFormat,
    source_name: str = "",
    selection: RecordSelection = "all",
    bonds: BondStrategy = "none",
    conformers: ConformerSelection = "all",
) -> MoleculeCollection:
    """Parse molecular text using an explicitly selected native reader."""

    if not isinstance(contents, str):
        raise TypeError("contents must be a string")
    if not isinstance(source_name, str):
        raise TypeError("source_name must be a string")
    _validate_options(format, selection, bonds, conformers)
    return _collection(
        _native_adapters._parse(contents, source_name, format, selection, bonds, conformers),
        source_name,
    )


def read(
    path: str | PathLike[str],
    *,
    format: InputFormat,
    selection: RecordSelection = "all",
    bonds: BondStrategy = "none",
    conformers: ConformerSelection = "all",
) -> MoleculeCollection:
    """Read UTF-8 molecular text using an explicitly selected native reader."""

    if not isinstance(path, (str, PathLike)):
        raise TypeError("path must be a string or path-like value")
    source_path = Path(path)
    return parse(
        source_path.read_text(encoding="utf-8"),
        format=format,
        source_name=str(source_path),
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


__all__ = [
    "InputFormat",
    "RecordSelection",
    "BondStrategy",
    "ConformerSelection",
    "parse",
    "read",
]
