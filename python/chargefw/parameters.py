"""Immutable parameter-set descriptors and explicit JSON loading."""

from __future__ import annotations

from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any

from ._chargefw import parameters as _native_parameters


@dataclass(frozen=True, slots=True)
class ParameterSetDescriptor:
    """Value-only metadata for an immutable parameter set."""

    id: str
    method_id: str
    name: str
    publication: str
    notes: str
    priority: int


class ParameterSet:
    """An owned immutable parameter set loaded from one JSON file."""

    __slots__ = ("_native", "_descriptor")

    def __init__(self, native: Any) -> None:
        descriptors = native._descriptors()
        if len(descriptors) != 1:
            raise RuntimeError("a Python ParameterSet must contain exactly one native parameter set")
        self._native = native
        self._descriptor = _descriptor(descriptors[0])

    @property
    def descriptor(self) -> ParameterSetDescriptor:
        return self._descriptor

    @property
    def id(self) -> str:
        return self._descriptor.id


def _descriptor(value: dict[str, Any]) -> ParameterSetDescriptor:
    return ParameterSetDescriptor(
        id=value["id"],
        method_id=value["method_id"],
        name=value["name"],
        publication=value["publication"],
        notes=value["notes"],
        priority=value["priority"],
    )


def _normalized_path(value: str | PathLike[str], name: str) -> str:
    if not isinstance(value, (str, PathLike)):
        raise TypeError(f"{name} must be a path")
    return str(Path(value))


def load_parameter_set(path: str | PathLike[str]) -> ParameterSet:
    """Load one immutable parameter set from a JSON file."""

    return ParameterSet(_native_parameters._load_parameter_set(_normalized_path(path, "path")))


def load_parameter_sets(directory: str | PathLike[str]) -> tuple[ParameterSet, ...]:
    """Load immutable JSON parameter sets from one directory."""

    path = _normalized_path(directory, "directory")
    return tuple(ParameterSet(native) for native in _native_parameters._load_parameter_sets(path))
