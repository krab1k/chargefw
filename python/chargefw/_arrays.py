"""Private NumPy normalization helpers."""

from __future__ import annotations

from collections.abc import Hashable, Iterable, Sequence
import numbers
from typing import Any

import numpy as np


def as_integer_array(
    value: Any, field: str, ndim: int, empty_1d_shape: tuple[int, ...] | None = None
) -> np.ndarray:
    try:
        array = np.asarray(value)
    except Exception as error:
        raise TypeError(f"{field} must be an integer array") from error

    if (
        array.ndim == 0
        and array.dtype.kind == "O"
        and isinstance(value, Iterable)
        and not isinstance(value, (str, bytes))
    ):
        try:
            array = np.asarray(tuple(value))
        except Exception as error:
            raise TypeError(f"{field} must be an integer array") from error

    if array.ndim != ndim:
        if empty_1d_shape is not None and array.ndim == 1 and array.size == 0:
            return np.empty(empty_1d_shape, dtype=np.int64)
        raise ValueError(f"{field} must have rank {ndim}, got {array.ndim}")
    if array.dtype.kind == "b":
        raise TypeError(f"{field} must contain integers, not booleans")
    if array.size == 0 and array.dtype.kind == "f":
        return np.empty(array.shape, dtype=np.int64)
    if array.dtype.kind == "O":
        if not all(
            isinstance(item, numbers.Integral) and not isinstance(item, (bool, np.bool_))
            for item in array.flat
        ):
            raise TypeError(f"{field} must contain integers")
        try:
            return np.array(array, dtype=np.int64, order="C", copy=True)
        except (OverflowError, TypeError, ValueError) as error:
            raise ValueError(
                f"{field} contains values outside the supported integer range"
            ) from error
    if array.dtype.kind not in "iu":
        raise TypeError(f"{field} must contain integers")

    if array.dtype.kind == "u" and np.any(array > np.iinfo(np.int64).max):
        raise ValueError(f"{field} contains values outside the supported integer range")
    try:
        return np.array(array, dtype=np.int64, order="C", copy=True)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{field} contains values outside the supported integer range") from error


def as_coordinates(value: Any, atom_count: int) -> tuple[np.ndarray, bool]:
    if value is None:
        return np.empty((0, atom_count, 3), dtype=np.float64), False
    try:
        array = np.asarray(value)
    except Exception as error:
        raise TypeError("coordinates must be a numeric array") from error
    if array.ndim == 2:
        if array.shape != (atom_count, 3):
            raise ValueError(f"coordinates must have shape ({atom_count}, 3)")
        one_conformer = True
        canonical_shape = (1, atom_count, 3)
    elif array.ndim == 3:
        if array.shape[1:] != (atom_count, 3):
            raise ValueError(f"coordinates must have shape (C, {atom_count}, 3)")
        one_conformer = False
        canonical_shape = array.shape
    else:
        raise ValueError("coordinates must have rank 2 or 3")
    if array.dtype.kind not in "iuf" or array.dtype.kind == "b":
        raise TypeError("coordinates must contain real numbers")
    try:
        canonical = np.array(array, dtype=np.float64, order="C", copy=True).reshape(canonical_shape)
    except (OverflowError, ValueError) as error:
        raise ValueError("coordinates contain values outside the supported range") from error
    if not np.all(np.isfinite(canonical)):
        raise ValueError("coordinates must contain only finite values")
    public = canonical.reshape((atom_count, 3)) if one_conformer else canonical
    return public, one_conformer


def as_names(value: Sequence[str] | None, count: int, field: str) -> tuple[str, ...]:
    if value is None:
        return ("",) * count
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{field} must be a sequence of strings")
    try:
        result = tuple(value)
    except TypeError as error:
        raise TypeError(f"{field} must be a sequence of strings") from error
    if len(result) != count:
        raise ValueError(f"{field} must contain {count} values")
    if not all(isinstance(item, str) for item in result):
        raise TypeError(f"{field} must contain only strings")
    return result


def as_ids(
    value: Iterable[Hashable] | None, count: int, field: str
) -> tuple[Hashable, ...]:
    if value is None:
        return tuple(range(count))
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{field} must be a sequence of hashable values")
    try:
        result = tuple(value)
    except TypeError as error:
        raise TypeError(f"{field} must be a sequence of hashable values") from error
    if len(result) != count:
        raise ValueError(f"{field} must contain {count} values")
    for item in result:
        try:
            hash(item)
        except TypeError as error:
            raise ValueError(f"{field} must contain only hashable values") from error
    return result


def readonly_copy(array: np.ndarray) -> np.ndarray:
    result = array.copy()
    result.setflags(write=False)
    return result
