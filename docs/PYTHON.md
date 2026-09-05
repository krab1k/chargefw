# Python package

The `chargefw` Python package provides owned molecule values, immutable NumPy arrays, method and parameter
catalogs, reusable assessment plans, calculation results, and Gemmi conversion over the native ChargeFW
engine.

The package is synchronous and in-process. Native molecule construction, assessment, and calculation
release the GIL, and independent calculations can run concurrently. Calculations can report structured
progress and support cooperative cancellation through a per-request observer.

## Installation status

ChargeFW currently builds Python wheels from the source tree but does not yet publish or qualify a binary
wheel matrix. Python 3.10 or newer is required. Runtime dependencies are NumPy 1.26 or newer and
Gemmi 0.7.4.

```bash
uv build --quiet --wheel
uv pip install --link-mode=copy --reinstall dist/chargefw-*.whl
```

The wheel contains the private native extension, required shared libraries, bundled parameter JSON,
Python modules, type declarations, and the `py.typed` marker. Parameter discovery uses package resources
and does not depend on the current directory or an environment variable.

For a direct CMake build, set `CHARGEFW_BUILD_PYTHON=ON`. The option defaults to `OFF` for ordinary native
builds.

## Basic calculation

```python
import chargefw

molecule = chargefw.Molecule(
    atomic_numbers=[8, 1, 1],
    formal_charges=[0, 0, 0],
    bonds=[[0, 1, 1], [0, 2, 1]],
    coordinates=[
        [0.0, 0.0, 0.0],
        [0.96, 0.0, 0.0],
        [-0.24, 0.93, 0.0],
    ],
    name="water",
)

result = chargefw.calculate(molecule, method="eem", execution="full")
print(result.assignments[0].values)
```

`calculate()` also accepts a `MoleculeCollection` or any iterable of `Molecule` values.

For a collection, `assignments_by_molecule` groups the source-ordered assignments for each input molecule.
Use it with the original collection for ordinary molecule-by-molecule iteration:

```python
molecules = chargefw.MoleculeCollection([first_molecule, second_molecule])
result = chargefw.calculate(molecules, method="qeq", execution="full")

for molecule, assignments in zip(molecules, result.assignments_by_molecule, strict=True):
    for assignment in assignments:
        print(molecule.name, assignment.conformer_index, assignment.values)
```

Geometry-dependent methods produce one assignment for each conformer; geometry-independent methods
produce one assignment with `conformer_index` set to `None`.

Retrieve one assignment directly by its source indices:

```python
second_conformer_charges = result.assignment(molecule=0, conformer=1).values
```

## Molecules and ownership

`Molecule` accepts:

- `atomic_numbers`: an integer vector of shape `(N,)`;
- `formal_charges`: an optional integer vector of shape `(N,)`, defaulting to zero;
- `bonds`: optional integer rows `(first_atom, second_atom, order)` of shape `(B, 3)`;
- `coordinates`: `(N, 3)` for one conformer or `(C, N, 3)` for multiple conformers;
- optional molecule, atom, and conformer names;
- a `SourceIdentity` or its `source_name`, `record_index`, and `record_id` fields; and
- optional hashable `atom_ids` for round-trip mapping.

Inputs are validated, normalized to C-contiguous `int64` or `float64` arrays, and copied. Public arrays
are read-only, so later mutation or destruction of the caller's arrays cannot change the molecule.
Coordinates are always exposed with shape `(C, N, 3)`; `None` and `(0, N, 3)` both mean no conformers.

Atomic numbers 1–100 are accepted because they are represented by the bundled periodic table; individual
methods and parameter sets may support a smaller subset. Only bond orders 1, 2, and 3 are accepted. Self
bonds, duplicate bonds, invalid indices, non-integral integer input, unsupported atomic numbers,
mismatched shapes are rejected. Non-finite coordinates are retained, but methods that require geometry
are inapplicable when coordinates are missing or non-finite.

`MoleculeCollection` is an immutable source-ordered sequence with an optional name. Omitted atom IDs
default to source-order integers.

## Methods and parameter sets

The package-level `chargefw.methods` and `chargefw.parameter_sets` values are immutable ordered mappings
with lookup by stable ID:

```python
eem = chargefw.methods["eem"]

for method_id, method in chargefw.methods.items():
    print(method_id, method.name, method.supports_cutoff)

for parameter_set in eem.parameter_sets.values():
    print(parameter_set.id, parameter_set.name)

iterations = chargefw.methods["peoe"].options["iters"]
print(iterations.default, iterations.minimum)
```

Method descriptors expose names, publication metadata, priority, coordinate requirements, reduced-mode
capabilities, options, and associated bundled parameter sets. The Python API does not expose native
method objects, parameter classifications, or custom parameter-catalog construction. The native
[parameter-set JSON reference](PARAMETERS.md) defines classifier behavior, including permissive matching.

## Assessment and calculation policy

`assess()` and automatic `calculate()` accept keyword-only policy arguments:

| Argument | Values and defaults |
| --- | --- |
| `method` | Method ID, `Method`, or `None` |
| `parameter_set` | Parameter-set ID, `ParameterSet`, or `None` |
| `options` | Option mapping for an explicitly selected method |
| `options_by_method` | Method-ID to option-mapping overrides for automatic selection |
| `parameter_matching` | `"strict"` (default) or `"permissive"` |
| `execution` | `"auto"` (default), `"full"`, `"cutoff"`, or `"cover"` |
| `radius` | Reduced radius; explicit cutoff/cover require at least `8.0` |
| `charge_correction` | `"uniform"`, `"none"`, or `None` |
| `cutoff_threshold` | Automatic full-to-cutoff threshold; default `20_000`, `None` is unlimited |
| `cover_threshold` | Automatic cutoff-to-cover threshold; default `80_000`, `None` is unlimited |
| `threads` | Non-negative oneTBB thread limit; omitted or `0` delegates to oneTBB |

Flat `options` require an explicit method. `options` and `options_by_method` cannot be combined. Automatic
execution accepts an optional radius override but not a charge-correction override. Explicit full
execution rejects radius and correction arguments; explicit cutoff and cover require a radius.

```python
assessment = chargefw.assess(
    molecule,
    method="eem",
    execution="auto",
)

for plan in assessment.plans:
    print(plan.method.id, plan.policy.mode, plan.warnings)
```

An `Assessment` contains priority-ordered reusable `plans`, structured `rejections`, a `default_plan`,
and assessment time in `seconds`. Each `Plan` retains the prepared native state and exposes its method,
parameter set, complete validated options, concrete execution policy, and warnings.

Execute an exact assessed plan without repeating preparation or parameter classification:

```python
plan = assessment.default_plan
if plan is not None:
    result = chargefw.calculate(molecule, plan, threads=1)
```

Plans are bound to the exact molecule objects and collection name used during assessment. Selection
arguments cannot be supplied with a plan. A plan is reusable after its `Assessment` is released and can
be used by independent concurrent calculations.

## Progress and cancellation

Subclass `CalculationObserver` and pass it to `calculate()` with the keyword-only `observer` argument:

```python
from threading import Event, Lock

import chargefw


class RecordingObserver(chargefw.CalculationObserver):
    def __init__(self) -> None:
        self.events: list[chargefw.CalculationProgress] = []
        self.cancellation = Event()
        self.lock = Lock()

    def on_progress(self, progress: chargefw.CalculationProgress) -> None:
        with self.lock:
            self.events.append(progress)

    def cancelled(self) -> bool:
        return self.cancellation.is_set()


observer = RecordingObserver()
result = chargefw.calculate(molecule, method="eem", observer=observer)
```

Progress phases are `"computation_started"`, `"computation_finished"`, `"target_started"`,
`"target_finished"`, and `"fragment_progress"`. A target is one molecule/conformer pair for a
geometry-dependent method and one molecule for a geometry-independent method. Fragment events are
throttled and are emitted only by cutoff and cover execution.

Each immutable `CalculationProgress` snapshot contains the phase, effective execution mode, method ID,
target indices and count, fragment completion and count, molecule and optional conformer indices, and
elapsed computation seconds. Fields unrelated to the current phase retain their default zero or `None`
value.

Callbacks may originate on oneTBB worker threads. The Python GIL is held during each callback, but
observers that update compound state should still use synchronization and must not rely on callback
arrival order during parallel execution. Callback exceptions are reported through Python's unraisable
exception hook and do not alter calculation control flow.

Returning true from `cancelled()` requests termination at the next cancellation checkpoint. The call
then raises `CalculationCancelledError`; its `result` has status `"cancelled"` and contains no partial
charge assignments. Assessment itself is not observed.

## Results and failures

`CalculationResult` exposes:

- `status`;
- immutable `assignments`;
- immutable `assignments_by_molecule`, aligned with the calculation input collection;
- requested policy in `requested`;
- detached executed provenance in `plan`;
- rejected alternatives and warnings;
- optional `failure_message`; and
- applicability and computation `timings`.

Each `ChargeAssignment` contains a newly owned, read-only, C-contiguous `float64` vector together with its
molecule index, optional conformer index, `SourceIdentity`, and atom IDs.
`numpy.asarray(assignment)` returns the charge vector.

Geometry-dependent methods return one assignment per conformer. Geometry-independent methods return one
assignment per molecule without conformer identity.

Invalid Python values and request arguments raise `TypeError` or `ValueError`. Failed calculations raise
a typed `ChargeFWError` subclass:

- `InvalidInputError`;
- `NoExecutablePlanError`;
- `NumericalFailureError`;
- `CalculationCancelledError`.

Each exception retains the complete result as `exception.result`, including status, rejections, failure
text, and timings.

## Molecular input

Serialized molecular formats are read through `chargefw.io`. Parsing functions accept text, while
corresponding `read_*()` functions accept string or path-like filesystem paths:

```python
import chargefw

molecules = chargefw.io.read_pdb(
    "structure.pdb",
    selection="all",
    bonds="hybrid",
    conformers="all",
)
```

All format functions return a `MoleculeCollection`, including formats that contain exactly one molecule:

- `parse_mol()` and `read_mol()`;
- `parse_sdf()` and `read_sdf()`;
- `parse_mol2()` and `read_mol2()`;
- `parse_molecule_json()` and `read_molecule_json()`;
- `parse_pdb()` and `read_pdb()`; and
- `parse_mmcif()` and `read_mmcif()`.

MOL, SDF, and MOL2 always import their format-defined single conformer. Molecule JSON accepts
`conformers="first"` or `"all"`. PDB and mmCIF additionally accept the structural selection and bond
options described below.

The base package depends on the upstream `gemmi` Python package. Importing `chargefw` does not import
Gemmi; import the object-conversion module explicitly:

```python
from chargefw.io import gemmi as chargefw_gemmi

structure_molecules = chargefw_gemmi.from_structure(structure, bonds="hybrid")
document_molecules = chargefw_gemmi.from_document(document)
```

The Gemmi integration serializes upstream objects through mmCIF text and then uses ChargeFW's compiled
reader. This keeps selection and bond behavior aligned with serialized input without sharing C++ objects
between extension modules.

Structural input defaults are `selection="all"`, `bonds="none"`, and `conformers="all"`. Bond choices are
`"none"`, `"explicit"`, `"templates"`, and `"hybrid"`; selection choices are `"all"`,
`"polymers-and-ligands"`, and `"polymers"`. Their language-independent import semantics are defined in
the [PDB and mmCIF format reference](FORMATS.md#pdb-and-mmcif-input).

The Python package currently has no molecular output API, RDKit or Biopython integration, chemistry
preparation API, or asynchronous job API. Current distribution and integration work is tracked in the root
[TODO](../TODO.md).
