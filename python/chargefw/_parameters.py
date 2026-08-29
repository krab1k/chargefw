"""Immutable bundled parameter-set descriptors and catalogs."""

from __future__ import annotations

from dataclasses import dataclass

from ._catalog import _Catalog
from ._payloads import ParameterSetDescriptorPayload


@dataclass(frozen=True, slots=True, repr=False)
class ParameterSet:
    """Metadata for an installed immutable parameter set."""

    id: str
    method: str
    name: str
    publication: str
    notes: str
    priority: int

    def __repr__(self) -> str:
        return f"{type(self).__name__}(id={self.id!r}, method={self.method!r}, name={self.name!r})"


class ParameterSetCatalog(_Catalog[ParameterSet]):
    """Installed parameter sets keyed by ID."""

    __slots__ = ()

    def __init__(self, values: tuple[ParameterSet, ...]) -> None:
        super().__init__(values, "parameter-set")

    def for_method(self, method: str) -> ParameterSetCatalog:
        """Return parameter sets associated with one method ID."""

        if not isinstance(method, str):
            raise TypeError("method must be a string")
        return ParameterSetCatalog(
            tuple(value for value in self.values() if value.method == method)
        )


def _parameter_set(value: ParameterSetDescriptorPayload) -> ParameterSet:
    return ParameterSet(
        id=value["id"],
        method=value["method_id"],
        name=value["name"],
        publication=value["publication"],
        notes=value["notes"],
        priority=value["priority"],
    )


ParameterSet.__module__ = "chargefw"
ParameterSetCatalog.__module__ = "chargefw"
