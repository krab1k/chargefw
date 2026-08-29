"""Method applicability and execution report values."""

from __future__ import annotations

from dataclasses import dataclass

from ._catalog import _Catalog
from ._chargefw import calculation as _native_calculation
from ._chargefw import methods as _native_methods
from ._payloads import MethodDescriptorPayload
from .parameters import ParameterSetCatalog

PrerequisiteIssueKind = _native_methods.PrerequisiteIssueKind
ExecutionAvailability = _native_methods.ExecutionAvailability
ExecutionIssueKind = _native_methods.ExecutionIssueKind
ExecutionMode = _native_calculation.ExecutionMode
MethodOptionType = _native_methods.MethodOptionType
for _enum in (
    PrerequisiteIssueKind,
    ExecutionAvailability,
    ExecutionIssueKind,
    MethodOptionType,
):
    _enum.__module__ = __name__


@dataclass(frozen=True, slots=True)
class PrerequisiteIssue:
    kind: PrerequisiteIssueKind
    message: str
    molecule_index: int | None = None
    atom_index: int | None = None
    bond_index: int | None = None
    conformer_index: int | None = None


@dataclass(frozen=True, slots=True)
class ExecutionIssue:
    kind: ExecutionIssueKind
    message: str
    molecule_index: int | None = None


@dataclass(frozen=True, slots=True)
class ExecutionAssessment:
    mode: ExecutionMode
    availability: ExecutionAvailability
    issues: tuple[ExecutionIssue, ...] = ()


@dataclass(frozen=True, slots=True)
class MethodOptionDescriptor:
    id: str
    description: str
    type: MethodOptionType
    default: bool | int | float | str
    choices: tuple[bool | int | float | str, ...]
    minimum: bool | int | float | str | None
    minimum_inclusive: bool
    maximum: bool | int | float | str | None
    maximum_inclusive: bool


class MethodOptionCatalog(_Catalog[MethodOptionDescriptor]):
    """Immutable ordered method options with lookup by ID."""

    __slots__ = ()

    def __init__(self, values: tuple[MethodOptionDescriptor, ...]) -> None:
        super().__init__(values, "method option")


@dataclass(frozen=True, slots=True)
class MethodDescriptor:
    id: str
    name: str
    full_name: str
    publication: str | None
    priority: int
    requires_coordinates: bool
    supports_cutoff: bool
    supports_cover: bool
    options: MethodOptionCatalog
    parameter_sets: ParameterSetCatalog


class MethodCatalog(_Catalog[MethodDescriptor]):
    """Immutable ordered built-in methods with lookup by ID."""

    __slots__ = ()

    def __init__(self, values: tuple[MethodDescriptor, ...]) -> None:
        super().__init__(values, "method")


def _method_descriptor(
    value: MethodDescriptorPayload, parameter_sets: ParameterSetCatalog
) -> MethodDescriptor:
    return MethodDescriptor(
        id=value["id"],
        name=value["name"],
        full_name=value["full_name"],
        publication=value["publication"],
        priority=value["priority"],
        requires_coordinates=value["requires_coordinates"],
        supports_cutoff=value["supports_cutoff"],
        supports_cover=value["supports_cover"],
        options=MethodOptionCatalog(
            tuple(
                MethodOptionDescriptor(
                    id=option["id"],
                    description=option["description"],
                    type=option["type"],
                    default=option["default"],
                    choices=tuple(option["choices"]),
                    minimum=option["minimum"],
                    minimum_inclusive=option["minimum_inclusive"],
                    maximum=option["maximum"],
                    maximum_inclusive=option["maximum_inclusive"],
                )
                for option in value["options"]
            )
        ),
        parameter_sets=parameter_sets.for_method(value["id"]),
    )


def _method_catalog(parameter_sets: ParameterSetCatalog) -> MethodCatalog:
    """Build calculator-bound method descriptors from the native registry."""

    return MethodCatalog(
        tuple(
            _method_descriptor(value, parameter_sets)
            for value in _native_methods._method_descriptors()
        )
    )
