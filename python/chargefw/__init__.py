"""ChargeFW Python API."""

from ._chargefw import version as _native_version

__version__ = _native_version()

__all__ = ["__version__"]
