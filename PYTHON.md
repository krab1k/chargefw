# ChargeFW Python Plan

This document owns the detailed Python scope, API contract, implementation plan, and delivery status.
Implemented cross-language architecture remains documented in [PROJECT.md](PROJECT.md), unfinished
release-level acceptance criteria in [TODO.md](TODO.md), and user instructions in
[README.md](README.md).

## Status

The Python package skeleton, optional nanobind target, owned array model, calculation vertical slice,
value-only method/parameter introspection, and required Gemmi integration are implemented. The bindings
are organized by the same domain boundaries as the native library. Calculation policy and report
vocabulary and Gemmi adapter policy use validated lowercase strings. No optional toolkit adapter is
implemented yet. The
native code already provides the foundations the binding should use rather than duplicate:

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
- inputs are copied into normalized C-contiguous `int64`/`float64` arrays backed by immutable storage;
  public array properties return those zero-copy read-only arrays, so caller mutation or lifetime cannot
  alter a molecule;
- validation rejects incorrect ranks/shapes, non-integral or overflowing integer values, unsupported
  atomic numbers and bond orders, self/duplicate/out-of-range bonds, non-finite coordinates, invalid
  names, negative record indices, and unhashable source IDs;
- omitted atom and conformer IDs default to source-order indices, and public coordinates always use
  canonical `(C, N, 3)` shape; `(0, N, 3)` represents no conformers while `(N, 3)` input is normalized
  to one conformer;
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

- `assess()` and `calculate()` accept validated keyword-only application policy, including flat
  options for an explicit method and method-scoped overrides for automatic selection;
- package-level `assess()` and `calculate()` load the bundled parameter catalog from package resources
  and accept a molecule, collection, or molecule iterable;
- `assess()` and `calculate()` use the native owned assessment facade, with
  assessment preparation and native calculation running while the GIL is released;
- `Assessment` exposes priority-ordered reusable plans, rejected scientific/policy alternatives, a
  default plan, and timing data while plans retain shared prepared native state;
- `CalculationResult` exposes status, source-mapped owned `float64` C-contiguous assignments, detached
  executed-plan provenance, execution issues, failure text, and timings; and
- no-plan, numerical-failure, cancellation, and invalid-request boundaries raise typed Python
  exceptions that retain their complete result;
- the native binding is split into `core`, `methods`, `parameters`, and `calculation` registration units,
  while public Python modules provide the corresponding high-level value API; and
- execution selection and charge correction use validated Python strings; effective modes, statuses,
  availability, option types, and issue kinds are normalized to documented lowercase strings.

Public plan, rejection, and provenance values use immutable tuples and read-only mappings. Plans retain
private shared prepared state so they remain independently executable after an `Assessment` is released;
prepared features, classifications, registry pointers, and other native internals are not exposed.

Validation added for this milestone:

```text
cmake --build build/python-skeleton-system --target chargefw_python --parallel 2
ctest --test-dir build/python-skeleton-system -R 'test_chargefw_python_(import|molecule|calculation)' --output-on-failure
uv build --quiet --wheel --out-dir /tmp/kilo/chargefw-m3-wheel
uv pip install --python /tmp/kilo/chargefw-m3-venv/bin/python --no-index --no-deps --reinstall \
  /tmp/kilo/chargefw-m3-wheel/chargefw-0.1-cp314-cp314-linux_x86_64.whl
```

The original binding structure was additionally validated with the focused Python suite, the complete
62-test `gcc-debug` native suite, a release wheel build, and an installed-wheel calculation. The public
policy surface has since been replaced by keyword-only arguments and lowercase strings.

The installed-wheel calculation smoke test ran from `/tmp` and succeeded. The offline environment did
not contain a NumPy wheel, so installation used `--no-deps` in a system-site-packages environment;
`pyproject.toml` now declares `numpy>=1.26` for normal online installation.

### Milestone 4 — bundled parameter introspection: complete

Implemented on 2026-08-27:

- package-level `methods`, each method's `options`, and package-level `parameter_sets` are immutable
  ordered mappings supporting standard keys, values, items, membership, `get()`, and lookup by stable
  ID; methods also expose their filtered bundled parameter sets;
- `ParameterSet` exposes immutable IDs, method IDs, names, publications, notes,
  and priorities; and
- the Python API deliberately exposes the installed bundled catalog only; arbitrary parameter tables,
  classifications, native method objects, and custom parameter-catalog construction remain private.

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

### Milestone 5 — Gemmi integration: complete

Implemented on 2026-08-28:

- the base Python distribution requires the tested upstream `gemmi==0.7.4` package without importing it
  from the top-level `chargefw` module;
- `chargefw.adapters.gemmi` accepts PDB/mmCIF text and files, `gemmi.Structure`, and
  `gemmi.cif.Document` values through serialization rather than cross-extension C++ object sharing;
- all serialized input is parsed by the existing native PDB/mmCIF readers, preserving native record
  selection, alternate-location handling, conformer selection, bond strategies, and failures;
- imported molecules preserve source record order and identity, atom names/IDs, model names/IDs,
  coordinates, formal charges, and supported bonds as ordinary owned Python values; and
- record selection (`all`, `polymers-and-ligands`, or `polymers`), bond strategy (`none`, `explicit`,
  `templates`, or `hybrid`), and conformer selection (`first` or `all`) use the same canonical lowercase
  strings as the CLI.

Validation performed:

```text
cmake --build --preset gcc-debug --parallel 2
ctest --preset gcc-debug --output-on-failure
cmake --build build/clang-debug --target chargefw_python --parallel 2
ctest --test-dir build/clang-debug -R 'test_chargefw_python_(import|gemmi|mypy)' --output-on-failure
uv build --quiet --wheel -Cbuild-dir=/tmp/kilo/chargefw-gemmi-build6 \
  --out-dir /tmp/kilo/chargefw-gemmi-wheel6
uv pip install --python /tmp/kilo/chargefw-gemmi-venv6/bin/python --link-mode=copy --reinstall \
  /tmp/kilo/chargefw-gemmi-wheel6/chargefw-0.1-cp314-cp314-linux_x86_64.whl
```

The complete 68-test GCC debug suite passed. Focused Clang debug tests, strict mypy, direct-component
installation, and an installed-wheel Gemmi/construction/calculation smoke test from `/tmp` also passed.
The wheel used CPython 3.14 on Linux x86-64 and nanobind 3.0.1. A missing nanobind STL string caster
include found during clean-wheel testing was corrected; ChargeFW remains compatible with both tested
nanobind 2.15 and 3.0 releases. Declaring and qualifying the release wheel matrix remains milestone 6.

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
  _methods.py              Private method metadata and mapping implementation
  _parameters.py           Private bundled parameter metadata implementation
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
submodules. The calculation boundary emits lowercase string vocabulary and accepts string policy values,
so native enum classes do not leak into the Python layer even privately. Whenever a binding registration
or its argument/return types change under `python/src/`, update the corresponding stub in the same change.
Keep these stubs implementation-focused; the public type contract belongs in the ordinary Python modules.
CMake copies the declarations into the build-tree package and installs them in wheels, while the import
smoke test verifies their presence.
The Clang ASan and UBSan presets build the extension and preload the matching Clang runtime for Python
CTest processes, because the host interpreter itself is not sanitizer-linked. The clang-tidy preset also
builds the extension so private binding sources are included in static analysis.

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

`coordinates=None` and shape `(0, N, 3)` both mean no conformers. The public `coordinates` property
always has shape `(C, N, 3)`, including `(1, N, 3)` for one conformer. Atom, molecule, conformer, and
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

result = chargefw.calculate(
    molecule,
    method="eem",
    execution="full",
)
charges = result.assignments[0].values
```

Calculation policy uses keyword-only arguments:

- optional method and parameter-set IDs or `Method`/`ParameterSet` values from the package mappings;
- flat `options` for an explicit method or advanced method-scoped `options_by_method` overrides;
- `parameter_matching="strict" | "permissive"`;
- `execution="auto" | "full" | "cutoff" | "cover"`, radius, and optional string charge correction;
- cutoff and cover atom thresholds, where `None` means unlimited; and
- calculation threads, where zero delegates to oneTBB and explicit limits must fit oneTBB's
  signed integer range.

Flat `options` require an explicit method. `options_by_method` preserves unambiguous overrides during
automatic selection, and the two forms cannot be combined.

The package loads its bundled parameter resources once and exposes immutable method and parameter-set
mappings. Read-only metadata exposes IDs, names, publications, priorities, option schemas, and
supported execution capabilities; parameter matching tables and native method objects remain private.
The initial Python API does not accept custom parameter catalogs.

The mappings preserve deterministic order and follow standard mapping behavior:

```python
eem = chargefw.methods["eem"]
iteration_option = chargefw.methods["peoe"].options["iters"]
eem_parameter_sets = chargefw.parameter_sets.for_method("eem")
same_parameter_sets = eem.parameter_sets

for method_id, method in chargefw.methods.items():
    print(method_id, method.name)
```

The following paths share one planning and execution implementation:

- `calculate(molecules, **policy)` assesses and executes the default plan;
- `assess(molecules, **policy)` returns priority-ordered runnable `plans`, `rejections`, and
  `default_plan`; and
- `calculate(molecules, plan)` executes that exact target-bound plan without repeating preparation,
  classification, or applicability checks. Plans are reusable and may be executed concurrently.

Explicit `full`, `cutoff`, or `cover` selection filters plans to that mode. Automatic selection exposes
all resource-policy-permitted modes and orders the preferred plan first. Plans excluded by automatic
resource limits remain visible as rejections; explicitly requested over-limit execution remains runnable
with a warning. A supplied plan cannot be combined with selection arguments and is rejected when passed
with a different molecule collection.

The extension releases the GIL during preparation and calculation. A later progress/cancellation
callback API may acquire the GIL for each callback, but it is not required for the first usable slice.
When added, callback exceptions must be captured, request cooperative cancellation, and be re-raised on
the calling Python thread; exceptions must never escape into oneTBB workers.

## Results, mapping, and failures

`CalculationResult` is an owned value. Reports, policies, issues, provenance, and timings are immutable
typed Python values rather than nested dictionaries. It exposes:

- lowercase status (`success`, `invalid_input_or_request`, `no_executable_plan`, `numerical_failure`, or
  `cancelled`);
- source-ordered assignments;
- rejected alternatives with structured issue kinds and indices;
- requested policy and detached executed-plan method/options/execution provenance;
- execution warnings and failure text; and
- applicability and computation timings.

Each assignment contains a newly owned, C-contiguous `float64` NumPy array plus molecule index,
optional conformer index, source record identity, atom IDs, and optional conformer ID.
Geometry-dependent methods return one assignment per molecule conformer. Geometry-independent methods
return one per molecule with no conformer identity. Returning Python-owned arrays avoids dangling views
when the native result or assessment is released. `ChargeAssignment` also implements NumPy's array
conversion protocol, so `np.asarray(assignment)` returns its immutable charge vector.

Construction and request-programming errors use normal Python exceptions: `TypeError` for incompatible
Python values and `ValueError` for invalid arrays, strings, options, or IDs. Scientific inapplicability,
numerical failure, and cancellation raise typed ChargeFW exceptions from `calculate()`; every exception
retains the complete result and diagnostics. Native exception text must retain molecule, conformer, and
method context.

## Adapter boundary

Adapters return ordinary `Molecule` or `MoleculeCollection` values. Calculation functions never detect
toolkit types or import optional toolkits. This keeps calculation policy in one place and lets adapters
evolve without changing the native facade.

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

Direct CMake installation places the Python component in the selected interpreter's platform package
directory under the install prefix (for example, `lib64/python3.14/site-packages/chargefw`), rather than
adding a nonstandard top-level directory. The extension and its required native runtime libraries remain
beside the package.

Development CMake presets pin `CHARGEFW_PYTHON_EXECUTABLE` to `/usr/bin/python3`. GCC and Clang debug and
release presets enable the bindings and their CTest suite. Wheel builds intentionally do not use those
presets: the PEP 517 frontend selects the target interpreter, which is necessary to build a wheel for each
supported CPython version.

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
needed; a clean build can select a new directory with `-Cbuild-dir=/tmp/chargefw-python-build`.

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

- `test_chargefw_python_mypy` runs strict mypy over the public package and binding tests, with the private
  extension stubs on `MYPYPATH`;
- accepted contiguous/non-contiguous array inputs and exact dtype/range/shape failures;
- input ownership after original arrays are mutated or destroyed;
- immutable molecule/collection ordering and source identity;
- C++ exception translation and contextual messages;
- bundled parameter mappings, option schemas, and deterministic selection;
- priority-ordered reusable plans, rejected alternatives, and no-plan results;
- full/cutoff/cover policy validation and detached executed-plan provenance;
- geometry-independent and all-conformer assignment cardinality;
- source molecule/atom/conformer mappings and Python-owned output arrays;
- repeated and concurrent functional calculations with the GIL released, including a
  worker-thread progress check during native work; and
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

### Performance conclusion

A one-time GCC release EEM comparison over matched native/Python inputs found no material binding latency:
at 2,000 atoms, native assessment/calculation medians were 1.25/62.6 ms and Python medians were
1.32/60.8 ms. The representative SDF corpus selected the same EEM parameter set in both facades; CLI phase
metrics were too coarsely rounded for comparison.

## Delivery sequence

Each milestone should leave one usable vertical slice and update this status section when complete.

1. **Complete — package and build skeleton** — add the optional CMake target, private extension, Python
   package, type marker, and focused import/version test; prove the native-only build is unchanged.
2. **Complete — owned array model** — implement molecule/collection construction, strict validation,
   source metadata, ownership tests, and conversion to native `core` values.
3. **Complete — calculation vertical slice** — bind keyword-only functional assess/calculate entry
   points, bundled parameter resources, NumPy assignments, reports, provenance, failures, and mappings.
4. **Complete — bundled catalog introspection** — expose value-only method/parameter metadata through
   immutable package mappings without exposing classification internals.
5. **Complete — Gemmi integration** — require the tested upstream Python package, implement serialized adapter
   entry points, and verify parity with native PDB/mmCIF import semantics.
6. **Wheel qualification** — build the declared initial matrix and run clean-environment relocation,
   shared-library, resource, and calculation smoke tests.
7. **RDKit integration** — add the lazy pure-Python converter and explicit attachment helper after the
   toolkit-neutral contract and wheel are stable.
8. **ACC III capability closure** — compare actual backend calls against the delivered API, add only
   missing explicit capabilities, and document intentional legacy exclusions before migration.

Biopython and additional distribution formats remain follow-on work driven by a concrete integration,
not prerequisites for the first Python release.
