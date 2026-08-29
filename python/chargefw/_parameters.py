"""Immutable bundled parameter-set descriptors and catalogs."""

from __future__ import annotations

from dataclasses import dataclass

from ._catalog import _Catalog
from ._payloads import ParameterSetDescriptorPayload


@dataclass(frozen=True, slots=True)
class ParameterSetDescriptor:
    """Value-only metadata for an immutable parameter set."""

    id: str
    method_id: str
    name: str
    publication: str
    notes: str
    priority: int


class ParameterSetCatalog(_Catalog[ParameterSetDescriptor]):
    """Immutable ordered parameter-set descriptors with ID and method lookup."""

    __slots__ = ()

    def __init__(self, values: tuple[ParameterSetDescriptor, ...]) -> None:
        super().__init__(values, "parameter-set")

    def for_method(self, method: str) -> ParameterSetCatalog:
        """Return parameter sets associated with one method ID."""

        if not isinstance(method, str):
            raise TypeError("method must be a string")
        return ParameterSetCatalog(tuple(value for value in self if value.method_id == method))


def _descriptor(value: ParameterSetDescriptorPayload) -> ParameterSetDescriptor:
    return ParameterSetDescriptor(
        id=value["id"],
        method_id=value["method_id"],
        name=value["name"],
        publication=value["publication"],
        notes=value["notes"],
        priority=value["priority"],
    )


ParameterSetDescriptor.__module__ = "chargefw"
ParameterSetCatalog.__module__ = "chargefw"
