"""Serialized integration with the upstream Gemmi Python package."""

from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import Literal, TypeAlias

import gemmi as _gemmi

from .._chargefw import adapters as _native_adapters
from ..core import Molecule, MoleculeCollection

_RecordSelection: TypeAlias = Literal["all", "polymers-and-ligands", "polymers"]
_BondStrategy: TypeAlias = Literal["none", "explicit", "templates", "hybrid"]
_ConformerSelection: TypeAlias = Literal["first", "all"]


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
        conformer_ids=payload["conformer_names"],
    )


def _path(value: str | PathLike[str]) -> Path:
    if not isinstance(value, (str, PathLike)):
        raise TypeError("path must be a string or path-like value")
    return Path(value)


def read_pdb_string(
    contents: str,
    *,
    source: str = "",
    selection: _RecordSelection = "all",
    bonds: _BondStrategy = "none",
    conformers: _ConformerSelection = "all",
) -> Molecule:
    """Read one molecule from PDB text using the native Gemmi adapter."""

    records = _native_adapters._read_pdb(contents, source, selection, bonds, conformers)
    if len(records) != 1:
        raise RuntimeError("PDB input did not produce exactly one molecule")
    return _molecule(records[0])


def read_pdb(
    path: str | PathLike[str],
    *,
    selection: _RecordSelection = "all",
    bonds: _BondStrategy = "none",
    conformers: _ConformerSelection = "all",
) -> Molecule:
    """Read one molecule from a UTF-8 PDB file."""

    source_path = _path(path)
    return read_pdb_string(
        source_path.read_text(encoding="utf-8"),
        source=str(source_path),
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


def read_mmcif_string(
    contents: str,
    *,
    source: str = "",
    selection: _RecordSelection = "all",
    bonds: _BondStrategy = "none",
    conformers: _ConformerSelection = "all",
) -> MoleculeCollection:
    """Read source-ordered molecules from mmCIF text using the native Gemmi adapter."""

    records = _native_adapters._read_mmcif(contents, source, selection, bonds, conformers)
    return MoleculeCollection((_molecule(record) for record in records), name=source)


def read_mmcif(
    path: str | PathLike[str],
    *,
    selection: _RecordSelection = "all",
    bonds: _BondStrategy = "none",
    conformers: _ConformerSelection = "all",
) -> MoleculeCollection:
    """Read source-ordered molecules from a UTF-8 mmCIF file."""

    source_path = _path(path)
    return read_mmcif_string(
        source_path.read_text(encoding="utf-8"),
        source=str(source_path),
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


def from_structure(
    structure: _gemmi.Structure,
    *,
    source: str = "",
    selection: _RecordSelection = "all",
    bonds: _BondStrategy = "none",
    conformers: _ConformerSelection = "all",
) -> Molecule:
    """Convert an upstream ``gemmi.Structure`` through mmCIF serialization."""

    if not isinstance(structure, _gemmi.Structure):
        raise TypeError("structure must be a gemmi.Structure")
    collection = read_mmcif_string(
        structure.make_mmcif_document().as_string(),
        source=source,
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )
    if len(collection) != 1:
        raise RuntimeError("Gemmi structure did not produce exactly one molecule")
    return collection[0]


def from_document(
    document: _gemmi.cif.Document,
    *,
    source: str = "",
    selection: _RecordSelection = "all",
    bonds: _BondStrategy = "none",
    conformers: _ConformerSelection = "all",
) -> MoleculeCollection:
    """Convert an upstream ``gemmi.cif.Document`` through mmCIF serialization."""

    if not isinstance(document, _gemmi.cif.Document):
        raise TypeError("document must be a gemmi.cif.Document")
    return read_mmcif_string(
        document.as_string(),
        source=source,
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


__all__ = [
    "read_pdb_string",
    "read_pdb",
    "read_mmcif_string",
    "read_mmcif",
    "from_structure",
    "from_document",
]
