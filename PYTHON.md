# ChargeFW Python Plan

This document owns the detailed Python scope, API contract, implementation plan, and delivery status.
Implemented cross-language architecture remains documented in [PROJECT.md](PROJECT.md), unfinished
release-level acceptance criteria in [TODO.md](TODO.md), and user instructions in
[README.md](README.md).

## Status

The Python package skeleton, optional nanobind target, owned array model, calculation vertical slice,
and value-only catalog introspection are implemented. The bindings are organized by the same domain
boundaries as the native library, and enum-valued policy and report fields use nanobind-backed Python
enums. No toolkit adapter is implemented yet. The native code already provides the foundations the
binding should use rather than duplicate:

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

The wheel smoke test currently covers CPython 3.14 on Linux x86-64 only. Gemmi metadata and full
clean-install qualification remain for later milestones.

### Milestone 2 — owned array model: complete

Implemented on 2026-08-27:

- `chargefw.Molecule` is an immutable Python value object accepting atomic numbers, formal charges,
  indexed bonds, one or more coordinate conformers, names, source identity, and source atom/conformer IDs;
- inputs are copied into normalized C-contiguous `int64`/`float64` arrays, with read-only copies returned
  by public array properties so caller mutation or lifetime cannot alter a molecule;
- validation rejects incorrect ranks/shapes, non-integral or overflowing integer values, unsupported
  atomic numbers and bond orders, self/duplicate/out-of-range bonds, non-finite coordinates, invalid
  names, negative record indices, and unhashable source IDs;
- omitted atom and conformer IDs default to source-order indices, and `(0, N, 3)` coordinates represent no
  conformers while `(N, 3)` represents one conformer;
- `MoleculeCollection` is an immutable source-ordered sequence with a collection name; and
- the Python owner converts normalized data once to private nanobind-backed native `core::Molecule`
  values; an owned native `core::MoleculeCollection` is materialized only for assessment, without adding
  Python metadata to the toolkit-neutral core.

Validation performed:

```text
cmake --build build/python-skeleton-system --target chargefw_python --parallel 2
ctest --test-dir build/python-skeleton-system -R 'test_chargefw_python_(import|molecule)' --output-on-failure
python3 -m compileall -q python/chargefw tests/python/test_molecule.py tests/python/test_import.py
clang-format --dry-run --Werror python/src/module.cpp
git diff --check
```

The complete CTest invocation in this build directory was not a valid full-suite result because many
pre-existing native test executables were configured but had not been built; both focused Python tests
passed. The next milestone can consume the private native objects through the owned calculation facade.

### Milestone 3 — calculation vertical slice: complete

Implemented on 2026-08-27:

- `CalculationOptions` validates application policy and preserves method-scoped option overrides;
- `Calculator` loads the bundled parameter catalog from package resources and accepts a molecule,
  collection, or molecule iterable;
- `Calculator.assess()` and `Calculator.calculate()` use the native owned assessment facade, with
  assessment preparation and native calculation running while the GIL is released;
- `Assessment` exposes copied applicability, execution-policy, warning, and timing data and enforces
  one-shot calculation ownership;
- `CalculationResult` exposes status, source-mapped owned `float64` C-contiguous assignments, reports,
  requested/effective provenance, execution issues, failure text, and timings; and
- no-plan, numerical-failure, cancellation, and invalid-request exception/result boundaries are wired
  through typed Python exceptions and result statuses;
- the native binding is split into `core`, `methods`, `parameters`, and `calculation` registration units,
  while public Python modules provide the corresponding high-level value API; and
- execution selection, effective execution mode, charge correction, status, availability, and issue
  kinds use native enum types rather than binding-local strings.

Public report and provenance values are immutable Python dataclasses containing tuples, read-only
mappings, and native enum values. They remain valid after native assessment ownership is consumed and
do not expose prepared features, classifications, registry pointers, or other native internals.

Validation added for this milestone:

```text
cmake --build build/python-skeleton-system --target chargefw_python --parallel 2
ctest --test-dir build/python-skeleton-system -R 'test_chargefw_python_(import|molecule|calculation)' --output-on-failure
uv build --quiet --wheel --out-dir /tmp/kilo/chargefw-m3-wheel
uv pip install --python /tmp/kilo/chargefw-m3-venv/bin/python --no-index --no-deps --reinstall \
  /tmp/kilo/chargefw-m3-wheel/chargefw-0.1-cp314-cp314-linux_x86_64.whl
```

The binding-structure and enum refinement was additionally validated with the focused Python suite,
the complete 62-test `gcc-debug` native suite, a release wheel build, and an installed-wheel calculation
using the domain modules and enum-valued policy API.

The installed-wheel calculation smoke test ran from `/tmp` and succeeded. The offline environment did
not contain a NumPy wheel, so installation used `--no-deps` in a system-site-packages environment;
`pyproject.toml` now declares `numpy>=1.26` for normal online installation.

### Milestone 4 — introspection and external parameters: complete

Implemented on 2026-08-27:

- `method_descriptors()` and `Calculator.methods` expose immutable value-only method IDs, names,
  publications, priorities, coordinate requirements, cutoff/cover capabilities, and complete option
  schemas, including native-backed option type enums;
- `ParameterSetDescriptor` exposes immutable parameter-set IDs, method IDs, names, publications, notes,
  and priorities;
- `load_parameter_set()` and `load_parameter_sets()` load immutable external JSON parameter values through
  the native parser without exposing classifications, parameter tables, or native method objects; and
- `Calculator(parameter_sets)` accepts an explicit non-empty, unique-ID sequence of loaded parameter sets
  which replaces the bundled catalog, while `Calculator.parameter_sets` exposes its value-only catalog.

Validation performed:

```text
cmake --build build/python-skeleton-system --target chargefw_python --parallel 2
ctest --test-dir build/python-skeleton-system -R 'test_chargefw_python_(import|molecule|calculation)' --output-on-failure
python3 -m compileall -q python/chargefw tests/python/test_import.py tests/python/test_calculation.py
clang-format --dry-run --Werror python/src/native_parameter_catalog.h python/src/parameters.cpp python/src/methods.cpp
git diff --check
```

The complete configured suite passed all ChargeFW tests, including the focused binding tests. Two existing
test-configuration failures remain unrelated to this milestone: Gemmi's unbuilt `cpptest` executable and
the custom-install-layout test's missing configured CLI11 source directory.

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
  __init__.py              Stable convenience re-exports and version
  core.py                  Molecule ownership and source identity
  calculation.py           Policy, assessment, execution, reports, and failures
  charges.py               Source-mapped charge assignments
  methods.py               Applicability and execution report values
  parameters.py            Parameter loading and descriptors (milestone 4)
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

### Private native stubs

`py.typed` is intentionally empty: it is the PEP 561 marker for the complete typed package. The checked-in
`python/chargefw/_chargefw/*.pyi` declarations describe the private nanobind extension and its native
submodules. Whenever a binding registration, its argument/return types, or a native enum changes under
`python/src/`, update the corresponding stub in the same change. Keep these stubs implementation-focused;
the public type contract belongs in the ordinary Python modules. CMake copies the declarations into the
build-tree package and installs them in wheels, while the import smoke test verifies their presence.

## Molecule contract

`chargefw.Molecule` is an immutable, owned Python-facing value. Its primary constructor accepts:

- one-dimensional integer `atomic_numbers` of length `N`;
- optional one-dimensional integer `formal_charges` of length `N`, defaulting to zero;
- optional indexed bonds as an integer array of shape `(B, 3)`, with rows
  `(first_atom_index, second_atom_index, order)` and supported orders 1, 2, or 3;
- optional coordinates of shape `(N, 3)` for one conformer or `(C, N, 3)` for `C` conformers;
- molecule, atom, and conformer names;
- a source identity containing source name, zero-based record index, and record ID; and
- optional atom and conformer IDs used only for round-trip mapping.

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
    chargefw.CalculationOptions(
        method="eem",
        execution=chargefw.ExecutionSelectionKind.FULL,
    ),
)
result.raise_for_status()
charges = result.assignments[0].values
```

`CalculationOptions` contains only application policy:

- optional method and parameter-set IDs;
- method-scoped option overrides as
  `Mapping[str, Mapping[str, bool | int | float | str]]`;
- permissive parameter typing;
- execution selection as `ExecutionSelectionKind`, radius, and optional `ChargeCorrectionPolicy`;
- cutoff and cover atom thresholds, where `None` means unlimited; and
- maximum calculation threads, where zero delegates to oneTBB and explicit limits must fit oneTBB's
  signed integer range.

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

`CalculationResult` is an owned value. Reports, policies, issues, provenance, and timings are immutable
typed Python values rather than nested dictionaries. It exposes:

- status as `ExecutionStatus` (`SUCCESS`, `INVALID_INPUT_OR_REQUEST`, `NO_EXECUTABLE_PLAN`,
  `NUMERICAL_FAILURE`, or `CANCELLED`);
- source-ordered assignments;
- applicable and rejected candidate reports with structured issue kinds and indices;
- requested and effective method/options/execution provenance;
- execution warnings and failure text; and
- applicability and computation timings.

Each assignment contains a newly owned, C-contiguous `float64` NumPy array plus molecule index,
optional conformer index, source record identity, atom IDs, and optional conformer ID.
Geometry-dependent methods return one assignment per molecule conformer. Geometry-independent methods
return one per molecule with no conformer identity. Returning Python-owned arrays avoids dangling views
when the native result or assessment is released.

Construction and request-programming errors use normal Python exceptions: `TypeError` for incompatible
Python values, including strings supplied for enum-valued fields; `ValueError` for invalid
arrays/options/IDs; and `IndexError` only for explicit sequence indexing. Scientific inapplicability,
numerical failure, and cancellation remain result statuses so their reports are not lost.
`result.raise_for_status()` provides typed ChargeFW exceptions carrying the same result for callers that
prefer exception flow. Native exception text must retain molecule, conformer, and method context.

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

Development CMake presets pin `CHARGEFW_PYTHON_EXECUTABLE` to `/usr/bin/python3`, so every local GCC and
Clang configuration uses the system interpreter. Wheel builds intentionally do not use those presets:
the PEP 517 frontend selects the target interpreter, which is necessary to build a wheel for each supported
CPython version.

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
2. **Complete — owned array model** — implement molecule/collection construction, strict validation,
   source metadata, ownership tests, and conversion to native `core` values.
3. **Complete — calculation vertical slice** — bind `Calculator`, bundled parameter resources, options,
   assess/calculate, NumPy assignments, reports, provenance, failures, and all-conformer mappings.
4. **Complete — introspection and external parameters** — expose value-only method/parameter descriptors
   and explicit immutable parameter loading needed by ACC III without exposing classification internals.
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
