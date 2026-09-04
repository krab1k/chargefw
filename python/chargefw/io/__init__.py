"""Molecular file-format input backed by ChargeFW's native readers."""

from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import Literal, TypeAlias

from .._chargefw import adapters as _native_adapters
from ..core import Molecule, MoleculeCollection

RecordSelection: TypeAlias = Literal["all", "polymers-and-ligands", "polymers"]
BondStrategy: TypeAlias = Literal["none", "explicit", "templates", "hybrid"]
ConformerSelection: TypeAlias = Literal["first", "all"]


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


def _single(
    payloads: list[_native_adapters.MoleculePayload], format_name: str
) -> Molecule:
    if len(payloads) != 1:
        raise RuntimeError(f"{format_name} input did not produce exactly one molecule")
    return _molecule(payloads[0])


def _path(value: str | PathLike[str]) -> Path:
    if not isinstance(value, (str, PathLike)):
        raise TypeError("path must be a string or path-like value")
    return Path(value)


def parse_mol(contents: str, *, source_name: str = "") -> Molecule:
    """Parse one molecule from MOL text."""

    return _single(_native_adapters._parse_mol(contents, source_name), "MOL")


def read_mol(path: str | PathLike[str]) -> Molecule:
    """Read one molecule from a UTF-8 MOL file."""

    source_path = _path(path)
    return parse_mol(source_path.read_text(encoding="utf-8"), source_name=str(source_path))


def parse_sdf(contents: str, *, source_name: str = "") -> MoleculeCollection:
    """Parse source-ordered molecules from SDF text."""

    return _collection(_native_adapters._parse_sdf(contents, source_name), source_name)


def read_sdf(path: str | PathLike[str]) -> MoleculeCollection:
    """Read source-ordered molecules from a UTF-8 SDF file."""

    source_path = _path(path)
    return parse_sdf(source_path.read_text(encoding="utf-8"), source_name=str(source_path))


def parse_mol2(contents: str, *, source_name: str = "") -> MoleculeCollection:
    """Parse source-ordered molecules from Tripos MOL2 text."""

    return _collection(_native_adapters._parse_mol2(contents, source_name), source_name)


def read_mol2(path: str | PathLike[str]) -> MoleculeCollection:
    """Read source-ordered molecules from a UTF-8 Tripos MOL2 file."""

    source_path = _path(path)
    return parse_mol2(source_path.read_text(encoding="utf-8"), source_name=str(source_path))


def parse_molecule_json(
    contents: str,
    *,
    source_name: str = "",
    conformers: ConformerSelection = "all",
) -> MoleculeCollection:
    """Parse ChargeFW molecule JSON 1.0 text."""

    return _collection(
        _native_adapters._parse_molecule_json(contents, source_name, conformers), source_name
    )


def read_molecule_json(
    path: str | PathLike[str], *, conformers: ConformerSelection = "all"
) -> MoleculeCollection:
    """Read a UTF-8 ChargeFW molecule JSON 1.0 document."""

    source_path = _path(path)
    return parse_molecule_json(
        source_path.read_text(encoding="utf-8"),
        source_name=str(source_path),
        conformers=conformers,
    )


def parse_pdb(
    contents: str,
    *,
    source_name: str = "",
    selection: RecordSelection = "all",
    bonds: BondStrategy = "none",
    conformers: ConformerSelection = "all",
) -> Molecule:
    """Parse one molecule from PDB text using Gemmi."""

    return _single(
        _native_adapters._parse_pdb(contents, source_name, selection, bonds, conformers), "PDB"
    )


def read_pdb(
    path: str | PathLike[str],
    *,
    selection: RecordSelection = "all",
    bonds: BondStrategy = "none",
    conformers: ConformerSelection = "all",
) -> Molecule:
    """Read one molecule from a UTF-8 PDB file."""

    source_path = _path(path)
    return parse_pdb(
        source_path.read_text(encoding="utf-8"),
        source_name=str(source_path),
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


def parse_mmcif(
    contents: str,
    *,
    source_name: str = "",
    selection: RecordSelection = "all",
    bonds: BondStrategy = "none",
    conformers: ConformerSelection = "all",
) -> MoleculeCollection:
    """Parse source-ordered molecules from mmCIF text using Gemmi."""

    return _collection(
        _native_adapters._parse_mmcif(contents, source_name, selection, bonds, conformers),
        source_name,
    )


def read_mmcif(
    path: str | PathLike[str],
    *,
    selection: RecordSelection = "all",
    bonds: BondStrategy = "none",
    conformers: ConformerSelection = "all",
) -> MoleculeCollection:
    """Read source-ordered molecules from a UTF-8 mmCIF file."""

    source_path = _path(path)
    return parse_mmcif(
        source_path.read_text(encoding="utf-8"),
        source_name=str(source_path),
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


__all__ = [
    "RecordSelection",
    "BondStrategy",
    "ConformerSelection",
    "parse_mol",
    "read_mol",
    "parse_sdf",
    "read_sdf",
    "parse_mol2",
    "read_mol2",
    "parse_molecule_json",
    "read_molecule_json",
    "parse_pdb",
    "read_pdb",
    "parse_mmcif",
    "read_mmcif",
]
