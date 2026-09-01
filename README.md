# ChargeFW

ChargeFW is a C++23 framework for empirical partial atomic-charge calculation. It provides a
toolkit-neutral native library, a command-line application, and Python bindings over the same
assessment and calculation engine.

The current implementation includes 22 charge methods, bundled parameter sets, deterministic
applicability and execution planning, source-ordered result mapping, native molecular-file adapters,
and full or explicit reduced calculations. It is not yet the production backend for Atomic Charge
Calculator III.

ChargeFW expects a molecular graph, formal charges, and any coordinates required by the selected
method. It does not parse SMILES, add hydrogens, assign protonation states, perceive arbitrary bonds,
or generate coordinates.

## Quick start

Build and install the native library and CLI:

```bash
cmake --preset gcc-release -DCMAKE_INSTALL_PREFIX="$PWD/_install"
cmake --build --preset gcc-release
cmake --install build/gcc-release --strip
```

Calculate charges for an SDF file:

```bash
_install/bin/chargefw calculate molecule.sdf output
```

Build a Python wheel from the source tree:

```bash
uv build --quiet --wheel
uv pip install --link-mode=copy --reinstall dist/chargefw-*.whl
```

```python
import chargefw

molecule = chargefw.Molecule(
    atomic_numbers=[8, 1, 1],
    bonds=[[0, 1, 1], [0, 2, 1]],
    coordinates=[[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
)
result = chargefw.calculate(molecule, method="eem", execution="full")
charges = result.assignments[0].values
```

The Python package is pre-release and does not yet have a qualified binary wheel matrix.

## Documentation

- [Project design and scientific scope](docs/PROJECT.md)
- [Command-line interface](docs/CLI.md)
- [Native C++ library](docs/NATIVE.md)
- [Python package](docs/PYTHON.md)
- [Unfinished work](TODO.md)
- [Repository contribution rules](AGENTS.md)

ChargeFW is distributed under the MIT license.
