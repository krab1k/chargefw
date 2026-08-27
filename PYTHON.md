# ChargeFW Python Plan

This document owns the detailed Python scope, API contract, implementation plan, and delivery status.
Implemented cross-language architecture remains documented in [PROJECT.md](PROJECT.md), unfinished
release-level acceptance criteria in [TODO.md](TODO.md), and user instructions in
[README.md](README.md).

## Status

The Python package skeleton, optional nanobind target, and wheel build are implemented. No usable
calculation API or toolkit adapter is implemented yet. The native code already provides the foundations
the binding should use rather than duplicate:

- owned toolkit-neutral molecules and conformers;
- an owned assessment/calculation facade with deterministic selection and structured applicability;
- source-ordered, molecule/conformer-indexed charge assignments and effective calculation provenance;
- bundled parameter loading relative to an installed native library; and
- required native Gemmi support for PDB/mmCIF import and mmCIF output.

The binding must expose the owned facade. Prepared features, parameter classifications, native method
pointers, low-level `CalculationInput`, and non-owning spans are implementation details and are not part
of the Python API.

### Milestone 1 — package and build skeleton: complete

Implemented on 2026-08-27:

- `CHARGEFW_BUILD_PYTHON` defaults to `OFF`, and `CHARGEFW_PYTHON_EXECUTABLE` can select the interpreter
  used for the binding independently of other CMake dependencies;
- CMake finds nanobind from that interpreter and builds the private `chargefw._chargefw` extension;
- scikit-build-core provides a PEP 517 build through `pyproject.toml`, including `uv build` support;
- the package skeleton exposes the native version and ships `py.typed`;
- Python-enabled installs place the ChargeFW, Gemmi, and oneTBB runtime libraries beside the extension
  with an origin-relative runtime path; and
- the focused import/version CTest passes while a native-only configuration builds with Python disabled.

Validation performed:

```text
uv build --wheel --out-dir /tmp/kilo/chargefw-wheel
uv pip install --python /tmp/kilo/chargefw-wheel-venv/bin/python --no-index --reinstall \
  /tmp/kilo/chargefw-wheel/chargefw-0.1-cp314-cp314-linux_x86_64.whl
cmake --build build/python-skeleton-system --target chargefw_python --parallel 2
ctest --test-dir build/python-skeleton-system -R test_chargefw_python_import --output-on-failure
cmake --build build/native-only-skeleton --target chargefw_core --parallel 2
```

The wheel smoke test currently covers CPython 3.14 on Linux x86-64 only. Runtime NumPy/Gemmi metadata,
parameter resources, molecule bindings, calculation bindings, and clean-install calculation tests remain
for later milestones.

## Product scope

Python is a first-class ChargeFW language, not a thin CLI wrapper. The initial package must support the
same scientific selection and execution policy as the native facade while presenting ordinary Python
value objects and NumPy arrays.

The initial supported workflow is synchronous in-process calculation:

1. construct one or more toolkit-neutral molecules from arrays, or convert them through an adapter;
2. assess applicability and resolve a concrete execution plan;
3. calculate without repeating preparation or classification; and
4. obtain source-mapped NumPy charge arrays, provenance, issues, and timing data.

The initial package does not provide chemistry preparation, SMILES parsing, bond perception,
protonation, hydrogen addition/removal, coordinate generation, geometry optimization, or an async/job
API. Toolkit adapters translate existing toolkit objects; they must not silently perform those tasks.

## Public package layout

```text
chargefw/
  __init__.py              Stable high-level API and version
  _chargefw.*              Private nanobind extension
  _data/parameters/*.json  Bundled parameter resources
  adapters/
    gemmi.py               Required Gemmi integration
    rdkit.py               Future lazy optional integration
    biopython.py           Future lazy optional integration
  py.typed                  Typing marker
```

Only documented names re-exported from `chargefw` and `chargefw.adapters` are public. The extension is
private so native organization and binding implementation can change without becoming a Python
compatibility promise. Type annotations and generated/native stubs are shipped with the wheel.

## Molecule contract

`chargefw.Molecule` is an immutable, owned Python-facing value. Its primary constructor accepts:

- one-dimensional integer `atomic_numbers` of length `N`;
- optional one-dimensional integer `formal_charges` of length `N`, defaulting to zero;
- optional indexed bonds as an integer array of shape `(B, 3)`, with rows
  `(first_atom_index, second_atom_index, order)` and supported orders 1, 2, or 3;
- optional coordinates of shape `(N, 3)` for one conformer or `(C, N, 3)` for `C` conformers;
- molecule, atom, and conformer names;
- a source identity containing source name, zero-based record index, and record ID; and
- optional source atom and conformer IDs used only for round-trip mapping.

Atom and conformer IDs may be Python hashable values. They remain in the Python owner and are returned
unchanged; the native calculation sees only source order and names. When IDs are omitted, zero-based
indices are used. This makes an adapter's mapping explicit without adding toolkit objects or opaque
Python references to `core::Molecule`.

Inputs are normalized and copied into native-owned memory during construction. NumPy views supplied by
the caller therefore need not outlive the molecule and may be mutable. Validation must reject, rather
than truncate or reinterpret:

- non-integral or overflowing atomic numbers, formal charges, indices, and bond orders;
- unsupported atomic numbers or bond orders;
- out-of-range/self/duplicate bonds according to the native molecule contract;
- wrong array rank or cardinality; and
- non-finite coordinates.

`coordinates=None` and shape `(0, N, 3)` both mean no conformers. Atom, molecule, conformer, and
collection order are always source order. Arbitrary array-like input may be accepted by the pure-Python
layer, but the normalized public properties are NumPy arrays with documented dtypes and shapes.

`MoleculeCollection` is a small immutable sequence carrying a collection name and molecules. Public
calculation functions accept either one `Molecule`, a sequence of molecules, or a collection and
normalize them once. Toolkit-specific metadata and import diagnostics remain attached to the
Python-facing molecule/collection and are copied into the corresponding result records.

## Calculation API

The target high-level surface is deliberately small:

```python
import numpy as np
import chargefw

molecule = chargefw.Molecule(
    atomic_numbers=np.array([8, 1, 1]),
    formal_charges=np.zeros(3, dtype=np.int8),
    bonds=np.array([[0, 1, 1], [0, 2, 1]]),
    coordinates=np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]]),
    name="water",
)

calculator = chargefw.Calculator()
result = calculator.calculate(
    molecule,
    chargefw.CalculationOptions(method="eem", execution="full"),
)
result.raise_for_status()
charges = result.assignments[0].values
```

`CalculationOptions` contains only application policy:

- optional method and parameter-set IDs;
- method-scoped option overrides as
  `Mapping[str, Mapping[str, bool | int | float | str]]`;
- permissive parameter typing;
- execution selection (`"auto"`, `"full"`, `"cutoff"`, or `"cover"`), radius, and correction;
- cutoff and cover atom thresholds, where `None` means unlimited; and
- maximum calculation threads, where zero delegates to oneTBB.

Method options remain method-scoped even when a method was selected explicitly. The Python layer must
not add context-sensitive shorthand that could apply an option to a different method after automatic
selection.

`Calculator` owns an immutable parameter-set catalog. With no argument it loads package resources once.
It also accepts an explicit sequence returned by `load_parameter_set()` or
`load_parameter_sets()`. An explicit sequence replaces, rather than silently merges with, the bundled
catalog. Read-only method and parameter-set descriptors expose IDs, names, publications, priorities,
option schemas, and supported execution capabilities; parameter matching tables and native method
objects remain private.

The following two paths share one implementation:

- `Calculator.calculate(molecules, options)` is the ordinary convenience path; and
- `Calculator.assess(molecules, options)` returns a one-shot `Assessment` with a value-only report,
  selected policy, and warnings. `Assessment.calculate()` consumes its native executable state without
  repeating work. The report remains readable after consumption; a second calculation raises
  `RuntimeError`.

The extension releases the GIL during preparation and calculation. A later progress/cancellation
callback API may acquire the GIL for each callback, but it is not required for the first usable slice.
When added, callback exceptions must be captured, request cooperative cancellation, and be re-raised on
the calling Python thread; exceptions must never escape into oneTBB workers.

## Results, mapping, and failures

`CalculationResult` is an owned value. It exposes:

- status (`success`, `invalid_input_or_request`, `no_executable_plan`, `numerical_failure`, or
  `cancelled`);
- source-ordered assignments;
- applicable and rejected candidate reports with structured issue kinds and indices;
- requested and effective method/options/execution provenance;
- execution warnings and failure text; and
- applicability and computation timings.

Each assignment contains a newly owned, C-contiguous `float64` NumPy array plus molecule index,
optional conformer index, source record identity, source atom IDs, and optional source conformer ID.
Geometry-dependent methods return one assignment per molecule conformer. Geometry-independent methods
return one per molecule with no conformer identity. Returning Python-owned arrays avoids dangling views
when the native result or assessment is released.

Construction and request-programming errors use normal Python exceptions: `TypeError` for incompatible
Python values, `ValueError` for invalid arrays/options/IDs, and `IndexError` only for explicit sequence
indexing. Scientific inapplicability, numerical failure, and cancellation remain result statuses so
their reports are not lost. `result.raise_for_status()` provides typed ChargeFW exceptions carrying the
same result for callers that prefer exception flow. Native exception text must retain molecule,
conformer, and method context.

## Adapter boundary

Adapters return ordinary `Molecule` or `MoleculeCollection` values. `Calculator` never detects toolkit
types or imports optional toolkits. This keeps calculation policy in one place and lets adapters evolve
without changing the native facade.

### Gemmi

Gemmi is different from optional chemistry toolkits: it is already a required public native dependency,
and Python builds must provide tested Gemmi Python integration. The base Python distribution therefore
depends on the upstream `gemmi` Python package at the native version tested by ChargeFW and provides
`chargefw.adapters.gemmi` functions for at least:

- PDB/mmCIF text and file input;
- `gemmi.Structure` input; and
- `gemmi.cif.Document` input.

The first implementation should cross from the upstream Gemmi module to ChargeFW through PDB/mmCIF
serialization and the existing native adapter, not by sharing C++ objects between independently built
extension modules. That boundary is slightly more expensive but avoids nanobind-domain, C++ ABI, symbol,
and duplicate-Gemmi ownership hazards. It also ensures Python and CLI structural selection, alternate
location, conformer, and bond-strategy semantics stay aligned. Profile before considering direct native
object interop.

Pin the tested Gemmi release exactly in the first wheel. Loosen the range only after adapter tests cover
multiple upstream releases. Source-tree binding builds remain optional and do not build or install a
second top-level `gemmi` module; Python adapter tests require the declared upstream package.

### RDKit

`chargefw.adapters.rdkit` is a pure-Python, lazy optional adapter. Importing `chargefw` must not import or
require RDKit. It converts selected conformers from `rdkit.Chem.Mol`, preserving atom indices, formal
charges, supported bond orders, names, and conformer IDs. Aromatic or otherwise unsupported bonds are
rejected unless an explicit, documented conversion policy is supplied. It performs no implicit
sanitization, hydrogen changes, protonation, embedding, or optimization.

Charge attachment is a separate helper. It writes only the assignment selected by the caller, records
method/provenance properties under a ChargeFW-owned namespace, and does not overwrite existing atom or
molecule properties unless `overwrite=True`.

### Biopython

The package layout reserves `chargefw.adapters.biopython`, but it is not an initial deliverable.
`Bio.PDB` objects do not reliably provide the complete formal-charge and bond graph required by many
methods. A future adapter may preserve structure/model/atom identities and coordinates, but it must
require explicit bond/formal-charge inputs when absent and must not infer chemistry. Add it only with a
concrete supported workflow and fixtures.

## Build and distribution

`CHARGEFW_BUILD_PYTHON` is an optional CMake option, defaulting to `OFF` for both top-level and subproject
native builds. When enabled it:

1. finds the selected Python interpreter and development module;
2. finds nanobind from the Python build environment;
3. builds the private `chargefw._chargefw` extension linked to `chargefw::core`; and
4. enables Python tests when `CHARGEFW_BUILD_TESTS=ON`.

Enabling Python must not change the C++ API, CLI, native install contents, or dependency behavior for a
normal build. Do not fetch Python or nanobind when the option is off.

Use `pyproject.toml` with scikit-build-core as the PEP 517 backend, nanobind as a build dependency, and
NumPy plus the tested upstream Gemmi package as runtime dependencies. The wheel contains the extension,
required native shared libraries, Python modules/stubs, and parameter JSON. Wheel repair tooling must
verify external shared-library policy and RPATHs; it must not rely on `LD_LIBRARY_PATH` or another
ChargeFW installation.

Python default parameter discovery uses `importlib.resources` to resolve `chargefw/_data/parameters`
and passes that explicit directory to the native loader. It does not use an environment variable or the
current working directory. Native-only installation continues to use the existing library-relative
discovery mechanism.

For local iteration, the scikit-build-core build directory is persistent and keyed by the wheel tag
(`build/python/{wheel_tag}`). CMake and FetchContent can therefore reuse the configured native core and
dependency build for repeated builds of the same interpreter/platform instead of rebuilding the unchanged
core in a temporary directory. A source change still causes the normal CMake dependency rebuild when
needed; a clean rebuild can use `-Ccmake.fresh=true` or remove the corresponding ignored build directory.

Typical low-output local commands are:

```text
uv build --quiet --wheel
uv pip install --link-mode=copy --reinstall dist/chargefw-*.whl
```

The `uv pip` hard-link warning is an installer warning, not a ChargeFW or wheel correctness problem. uv
normally hard-links packages from its cache into the target environment. If the cache and environment are
on different filesystems, hard-linking is impossible and uv safely falls back to copying. Use
`--link-mode=copy` for a deliberate warning-free copy, or set `UV_LINK_MODE=copy` for the shell/project
workflow. This affects installation speed and disk usage only, not the built wheel.

Start with a declared Linux x86-64 CPython matrix that the existing native toolchain can support; the
implementation milestone must record exact CPython and manylinux tags after a wheel-toolchain spike.
Do not claim macOS, Windows, musllinux, PyPy, or additional architectures until clean-install and native
dependency tests run on them.

## Validation

### Binding tests

- accepted contiguous/non-contiguous array inputs and exact dtype/range/shape failures;
- input ownership after original arrays are mutated or destroyed;
- immutable molecule/collection ordering and source identity;
- C++ exception translation and contextual messages;
- bundled and explicit parameter catalogs, option schemas, and deterministic selection;
- applicability/no-plan reports and one-shot assessment ownership;
- full/cutoff/cover policy validation and effective provenance;
- geometry-independent and all-conformer assignment cardinality;
- source molecule/atom/conformer mappings and Python-owned output arrays;
- repeated and concurrent calculations with the GIL released; and
- cancellation/progress callback safety when that API is introduced.

### Adapter tests

- Gemmi PDB/mmCIF text, `Structure`, and `cif.Document` equivalence;
- structural selection, alternate locations, model ordering, source IDs, and bond strategies;
- lazy failure with an actionable import error for absent optional toolkits;
- RDKit index/formal-charge/bond/conformer preservation and no implicit preparation; and
- charge attachment selection and overwrite protection.

### Package tests

Build the wheel, install it into a clean virtual environment, change to an unrelated working directory,
and verify that it can import `chargefw` and `gemmi`, discover bundled parameters without environment
variables, calculate a fixture, preserve mappings, and load required shared libraries. Also test an
editable install and a native build with `CHARGEFW_BUILD_PYTHON=OFF` and no Python development package.

## Delivery sequence

Each milestone should leave one usable vertical slice and update this status section when complete.

1. **Complete — package and build skeleton** — add the optional CMake target, private extension, Python
   package, type marker, and focused import/version test; prove the native-only build is unchanged.
2. **Owned array model** — implement molecule/collection construction, strict validation, source
   metadata, ownership tests, and conversion to native `core` values.
3. **Calculation vertical slice** — bind `Calculator`, bundled parameter resources, options,
   assess/calculate, NumPy assignments, reports, provenance, failures, and all-conformer mappings.
4. **Introspection and external parameters** — expose value-only method/parameter descriptors and
   explicit immutable parameter loading needed by ACC III without exposing classification internals.
5. **Gemmi integration** — require the tested upstream Python package, implement serialized adapter
   entry points, and verify parity with native PDB/mmCIF import semantics.
6. **Wheel qualification** — build the declared initial matrix and run clean-environment relocation,
   shared-library, resource, and calculation smoke tests.
7. **RDKit integration** — add the lazy pure-Python converter and explicit attachment helper after the
   toolkit-neutral contract and wheel are stable.
8. **ACC III capability closure** — compare actual backend calls against the delivered API, add only
   missing explicit capabilities, and document intentional legacy exclusions before migration.

Biopython and additional distribution formats remain follow-on work driven by a concrete integration,
not prerequisites for the first Python release.
