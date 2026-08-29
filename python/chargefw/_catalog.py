"""Private immutable ordered catalog support."""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import Generic, Protocol, TypeVar, overload


class _Identified(Protocol):
    @property
    def id(self) -> str: ...


T = TypeVar("T", bound=_Identified)


class _Catalog(Sequence[T], Generic[T]):
    """An immutable sequence with lookup by stable string ID."""

    __slots__ = ("_by_id", "_kind", "_values")

    _by_id: dict[str, T]
    _kind: str
    _values: tuple[T, ...]

    def __init__(self, values: Sequence[T], kind: str) -> None:
        ordered = tuple(values)
        by_id = {value.id: value for value in ordered}
        if len(by_id) != len(ordered):
            raise ValueError(f"{kind} IDs must be unique")
        object.__setattr__(self, "_values", ordered)
        object.__setattr__(self, "_by_id", by_id)
        object.__setattr__(self, "_kind", kind)

    def __setattr__(self, name: str, value: object) -> None:
        raise AttributeError(f"{type(self).__name__} is immutable")

    def __len__(self) -> int:
        return len(self._values)

    def __iter__(self) -> Iterator[T]:
        return iter(self._values)

    @overload
    def __getitem__(self, key: int) -> T: ...

    @overload
    def __getitem__(self, key: slice) -> tuple[T, ...]: ...

    @overload
    def __getitem__(self, key: str) -> T: ...

    def __getitem__(self, key: int | slice | str) -> T | tuple[T, ...]:
        if isinstance(key, str):
            try:
                return self._by_id[key]
            except KeyError:
                raise KeyError(f"unknown {self._kind} ID: {key!r}") from None
        return self._values[key]

    def __contains__(self, value: object) -> bool:
        if isinstance(value, str):
            return value in self._by_id
        return value in self._values

    def get(self, id: str, default: T | None = None) -> T | None:
        """Return the value with *id*, or *default* when it is absent."""

        if not isinstance(id, str):
            raise TypeError("catalog IDs must be strings")
        return self._by_id.get(id, default)

    def ids(self) -> tuple[str, ...]:
        """Return IDs in deterministic catalog order."""

        return tuple(value.id for value in self._values)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({list(self._values)!r})"
