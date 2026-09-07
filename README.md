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

## Python usage

### From a molecular file

```python
import chargefw

molecules = chargefw.io.read("molecule.sdf", format="sdf")
result = chargefw.calculate(
    molecules,
    method="qeq",
    parameter_set="QEq_original",
    execution="full",
)
print(result.assignments[0].values)
```

### From arrays

```python
import chargefw

molecule = chargefw.Molecule(
    atomic_numbers=[8, 1, 1],
    bonds=[[0, 1, 1], [0, 2, 1]],
    coordinates=[[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
)
result = chargefw.calculate(
    molecule,
    method="qeq",
    parameter_set="QEq_original",
    execution="full",
)
charges = result.assignments[0].values
```

The Python package is pre-release and does not yet have a qualified binary wheel matrix.

## Documentation

- **[Python API](docs/PYTHON.md)** - read molecules, calculate charges, and inspect methods and results.
- **[Python recipes](docs/recipes/README.md)** - run complete SDF, Gemmi, RDKit, and parameter-set
  workflows.
- **[Command-line interface](docs/CLI.md)** - calculate and inspect molecular files from the CLI.
- **[Molecular formats](docs/FORMATS.md)** - check supported input, output, parsing, and mapping
  semantics.
- **[Parameter sets](docs/PARAMETERS.md)** - understand parameter data and matching behavior.
- **[Native C++ API](docs/NATIVE.md)** - embed ChargeFW in a C++ application.
- **[Project design](docs/PROJECT.md)** - understand the architecture and scientific scope.
- **[Unfinished work](TODO.md)** - see current product gaps.
- **[Contribution guide](AGENTS.md)** - work on the ChargeFW repository.

ChargeFW is distributed under the MIT license.
