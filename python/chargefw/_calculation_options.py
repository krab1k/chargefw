"""Private calculation option normalization and validation."""

from __future__ import annotations

import numbers
from collections.abc import Mapping
from dataclasses import dataclass, field
from operator import index as as_index
from types import MappingProxyType
from typing import Any

import numpy as np

from ._chargefw import calculation as _native_calculation
from ._types import MethodOptionValue

_MAX_NATIVE_THREADS = int(np.iinfo(np.int32).max)


def _normalized_nonnegative_integer(value: Any, field_name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{field_name} must be an integer")
    try:
        result = as_index(value)
    except TypeError as error:
        raise TypeError(f"{field_name} must be an integer") from error
    if result < 0:
        raise ValueError(f"{field_name} must be non-negative")
    if result > np.iinfo(np.uintp).max:
        raise ValueError(f"{field_name} is outside the native size range")
    return result


def _normalized_method_options(
    values: Mapping[str, Mapping[str, MethodOptionValue]],
) -> dict[str, dict[str, MethodOptionValue]]:
    if not isinstance(values, Mapping):
        raise TypeError("method_options must be a mapping")
    result: dict[str, dict[str, MethodOptionValue]] = {}
    int_limits = np.iinfo(np.int32)
    for method_id, overrides in values.items():
        if not isinstance(method_id, str):
            raise TypeError("method option method IDs must be strings")
        if not isinstance(overrides, Mapping):
            raise TypeError(f"method_options[{method_id!r}] must be a mapping")
        normalized: dict[str, MethodOptionValue] = {}
        for option_id, value in overrides.items():
            if not isinstance(option_id, str):
                raise TypeError("method option IDs must be strings")
            if isinstance(value, (bool, np.bool_)):
                normalized[option_id] = bool(value)
            elif isinstance(value, numbers.Integral):
                integer = int(value)
                if integer < int_limits.min or integer > int_limits.max:
                    raise ValueError(f"method option {method_id}.{option_id} is outside int range")
                normalized[option_id] = integer
            elif isinstance(value, numbers.Real):
                floating = float(value)
                if not np.isfinite(floating):
                    raise ValueError(f"method option {method_id}.{option_id} must be finite")
                normalized[option_id] = floating
            elif isinstance(value, str):
                normalized[option_id] = value
            else:
                raise TypeError(
                    f"method option {method_id}.{option_id} must be bool, int, float, or str"
                )
        result[method_id] = normalized
    return result


def _frozen_method_options(
    values: Mapping[str, Mapping[str, MethodOptionValue]],
) -> Mapping[str, Mapping[str, MethodOptionValue]]:
    normalized = _normalized_method_options(values)
    return MappingProxyType(
        {
            method_id: MappingProxyType(overrides)
            for method_id, overrides in normalized.items()
        }
    )


@dataclass(frozen=True, slots=True)
class CalculationOptions:
    """Application policy for one synchronous calculation request."""

    method: str | None = None
    parameter_set_id: str | None = None
    method_options: Mapping[str, Mapping[str, MethodOptionValue]] = field(default_factory=dict)
    permissive_types: bool = False
    execution: _native_calculation.ExecutionSelectionKind = (
        _native_calculation.ExecutionSelectionKind.AUTOMATIC
    )
    radius: float | None = None
    charge_correction: _native_calculation.ChargeCorrectionPolicy | None = None
    cutoff_atom_threshold: int | None = 20_000
    cover_atom_threshold: int | None = 80_000
    max_threads: int = 0

    def __post_init__(self) -> None:
        for field_name in ("method", "parameter_set_id"):
            value = getattr(self, field_name)
            if value is not None and not isinstance(value, str):
                raise TypeError(f"{field_name} must be a string or None")
        if not isinstance(self.permissive_types, (bool, np.bool_)):
            raise TypeError("permissive_types must be a boolean")
        if not isinstance(self.execution, _native_calculation.ExecutionSelectionKind):
            raise TypeError("execution must be an ExecutionSelectionKind")
        if self.charge_correction is not None and not isinstance(
            self.charge_correction, _native_calculation.ChargeCorrectionPolicy
        ):
            raise TypeError("charge_correction must be a ChargeCorrectionPolicy or None")
        if (
            self.execution is _native_calculation.ExecutionSelectionKind.AUTOMATIC
            and self.charge_correction is not None
        ):
            raise ValueError("automatic execution does not accept a charge correction")
        if self.execution is _native_calculation.ExecutionSelectionKind.FULL:
            if self.radius is not None:
                raise ValueError("full execution does not accept a radius")
            if self.charge_correction is not None:
                raise ValueError("full execution does not accept a charge correction")
        if self.execution in (
            _native_calculation.ExecutionSelectionKind.CUTOFF,
            _native_calculation.ExecutionSelectionKind.COVER,
        ) and self.radius is None:
            raise ValueError(f"{self.execution.name.lower()} execution requires a radius")
        if self.radius is not None:
            if not isinstance(self.radius, numbers.Real) or isinstance(
                self.radius, (bool, np.bool_)
            ):
                raise TypeError("radius must be a real number or None")
            if not np.isfinite(float(self.radius)) or float(self.radius) < 8.0:
                raise ValueError("radius must be finite and at least 8.0")
        cutoff_threshold = (
            None
            if self.cutoff_atom_threshold is None
            else _normalized_nonnegative_integer(
                self.cutoff_atom_threshold, "cutoff_atom_threshold"
            )
        )
        cover_threshold = (
            None
            if self.cover_atom_threshold is None
            else _normalized_nonnegative_integer(
                self.cover_atom_threshold, "cover_atom_threshold"
            )
        )
        if cover_threshold is not None and cutoff_threshold is None:
            raise ValueError("cover_atom_threshold requires a finite cutoff_atom_threshold")
        if (
            cover_threshold is not None
            and cutoff_threshold is not None
            and cover_threshold < cutoff_threshold
        ):
            raise ValueError("cover_atom_threshold must not be smaller than cutoff_atom_threshold")
        max_threads = _normalized_nonnegative_integer(self.max_threads, "max_threads")
        if max_threads > _MAX_NATIVE_THREADS:
            raise ValueError("max_threads exceeds oneTBB's supported integer range")

        object.__setattr__(self, "permissive_types", bool(self.permissive_types))
        object.__setattr__(self, "method_options", _frozen_method_options(self.method_options))
        if self.radius is not None:
            object.__setattr__(self, "radius", float(self.radius))
        object.__setattr__(self, "cutoff_atom_threshold", cutoff_threshold)
        object.__setattr__(self, "cover_atom_threshold", cover_threshold)
        object.__setattr__(self, "max_threads", max_threads)
