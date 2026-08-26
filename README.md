# ChargeFW

ChargeFW is a C++23, library-first framework for empirical partial atomic-charge calculation and a
modern successor to ChargeFW2, the engine used by Atomic Charge Calculator III (ACC III).

The repository includes a toolkit-neutral C++ library, 22 built-in methods, parameter loading and
classification, applicability/execution planning, full and parallel reduced calculation, molecular-file
adapters, and a focused CLI. It is not yet a production ACC III backend, Python package, or general
SMILES/chemistry-preparation tool.

A bounded ChargeFW2/publication audit found numerical parity for the supported finite-molecule variants.
Periodic, Ewald, and other unsupported publication branches remain out of scope; reduced modes are explicit
approximations without a general accuracy claim.

## Documentation

| Document | Purpose |
| --- | --- |
| [PROJECT.md](PROJECT.md) | Implemented architecture, capabilities, limits, and product direction |
| [TODO.md](TODO.md) | Unfinished deliverables and acceptance criteria |
| [AGENTS.md](AGENTS.md) | Repository boundaries and implementation rules |

## Requirements

- CMake 3.28 or newer
- Ninja
- GCC or Clang with C++23 support
- Internet access on first configure unless FetchContent sources are already cached or supplied

CMake fetches pinned CLI11 2.7.2, nlohmann/json 3.12.0, Eigen 5.0.1, nanoflann 1.12.1, oneTBB 2023.1.0,
and Gemmi 0.7.4 into the build tree by default. Set `CHARGEFW_USE_SYSTEM_DEPENDENCIES=ON` to prefer
compatible installed packages and fetch only missing dependencies. When `CHARGEFW_BUILD_TESTS` is
enabled, CMake similarly searches for Snitch 1.3.2 and fetches it only if unavailable.

## Build and test

Run from the repository root:

```bash
cmake --preset gcc-debug
cmake --build --preset gcc-debug
ctest --preset gcc-debug
```

Run one focused test after building:

```bash
ctest --test-dir build/gcc-debug -R test_reduced_execution --output-on-failure
```

AddressSanitizer configuration:

```bash
cmake --preset clang-asan
cmake --build --preset clang-asan --parallel 1
ctest --preset clang-asan
```

Sanitized Clang builds can require substantial memory for template-heavy translation units; keep the
build serial (or use a deliberately chosen low `--parallel` value) to avoid exhausting system memory.
`clang-asan` enables only AddressSanitizer. Use the separate `clang-ubsan` preset for undefined-behavior
checks.

The tracked preset palette is `gcc-debug`, `clang-debug`, `gcc-release`, `clang-release`, `clang-asan`,
`clang-ubsan`, and `clang-tidy`. Useful options include `CHARGEFW_BUILD_TESTS`,
`CHARGEFW_BUILD_CLI`, `CHARGEFW_ENABLE_CCACHE`, `CHARGEFW_ENABLE_ASAN`,
`CHARGEFW_ENABLE_UBSAN`, `CHARGEFW_ENABLE_CLANG_TIDY`, and
`CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS` (enabled by default; disable for portable builds).

## CLI quick start

Install ChargeFW before running the CLI; bundled parameter JSON is resolved relative to the installed
library:

```bash
cmake --install build/gcc-debug --prefix "$PWD/_install"
_install/bin/chargefw calculate \
  tests/fixtures/synthetic/sdf/water.sdf output
```

Commands:

```text
chargefw calculate [options] INPUT OUTPUT_DIRECTORY
chargefw inspect [input-options] INPUT
chargefw applicability [options] INPUT
chargefw methods
chargefw parameters [METHOD]
```

Use `chargefw COMMAND --help` for the complete option list.

### Calculation and applicability options

- `--method ID` and `--parameter-set ID` restrict selection; omitted IDs use deterministic priorities.
- `--permissive-types` enables permissive parameter classification.
- `--execution auto|full|cutoff|cover` selects execution (`auto` is default).
- `--radius ANGSTROM` is required for explicit cutoff/cover and must be at least 8 Å. With automatic
  reduced execution it overrides the 12 Å default.
- `--charge-correction uniform|none` applies only to explicit reduced execution; uniform is the
  reduced-execution default.
- `--cutoff-atom-threshold COUNT|unlimited` changes the automatic full-to-cutoff threshold (default
  20,000 atoms for expensive full methods).
- `--cover-atom-threshold COUNT|unlimited` changes the automatic cutoff-to-cover threshold (default
  80,000 atoms).
- `--threads COUNT` limits calculation concurrency; `0` uses the oneTBB default.
- `--progress` displays a live calculation progress bar on standard error.
- `--method-option METHOD.OPTION=VALUE` supplies a repeatable, method-scoped option override. For
  automatic method selection, include the method ID; `chargefw methods` lists option schemas,
  defaults, and choices.

Automatic selection uses full execution when it is not discouraged. For an expensive candidate above
the cutoff threshold, it uses supported cutoff at 12 Å. Above the cover threshold, it uses supported
cover. Explicit execution overrides these thresholds; explicit full or cutoff records any applicable
threshold warning.

Cutoff and cover are available for ABEEM, EEM, EQeq, EQeq+C, QEq, SFKEEM, SQE, SQE+q0, and SQE+qp.
Cover retains source-ordered charges within 3 Å of each solved pivot halo; it remains a new
approximation without a general accuracy envelope.

### Input and output

Supported input extensions:

```text
.mol  .sdf  .mol2  .json  .pdb  .cif  .mmcif
```

Input format is selected by extension. The CLI rejects the collection on the first malformed record.
For PDB/mmCIF input:

- `--structural-selection all|polymers-and-ligands|polymers` (default `all`)
- `--structural-bonds none|explicit|templates|hybrid` (CLI default `hybrid`)

These structural options are rejected for other formats.

For all input formats, `--conformers first|all` selects whether the first conformer/model or every
conformer/model is read (default `all`). The option is ignored when a molecule has no conformers.
With `first` on PDB/mmCIF input, structural output may retain the complete source structure, but
ChargeFW writes charges only for the first calculated conformer and maps them to the corresponding
source atom IDs. JSON output records the selection in
`calculation_provenance.requested.input.conformers`.

Successful nonstructural calculations write:

```text
OUTPUT_DIRECTORY/<input-stem>.chargefw.{json,sdf,mol2,cif}
```

PDB/mmCIF input writes JSON and mmCIF only. Same-format SDF/MOL2 output preserves source content where
possible; other molecular output is generated. JSON input with multiple conformers cannot currently be
written to SDF/MOL2 and is rejected; use `--conformers first` for those output formats.

Completed calculations exit `0`; unexpected internal failures exit `1`. Invalid input or requests exit
`2`, no executable plan exits `3`, numerical calculation failure exits `4`, and cancellation exits `5`.
After input import, every outcome writes `<input-stem>.chargefw.json` to the requested directory; only
successful calculations write molecular charge files.

`inspect` reports imported composition without loading parameters. `applicability` reports scientific
candidates and their full/cutoff/cover availability without calculation. `methods` and `parameters`
list installed capabilities.

### JSON result

Schema `1.0` contains invocation and source-record statuses, diagnostics, source identity, and, only
on success, source-ordered charge assignments and totals. `calculation_provenance` records requested
conformer selection, selection/classification/resource/structural policies, and the effective method,
parameter set, method options, execution mode, radius, correction, and warnings. Requested method
options contain explicit overrides; effective options contain the complete selected-method values.
JSON charge values are rounded to at most four decimal places, and assignment totals sum the serialized
values; internal calculations retain native precision. Diagnostic messages use one-based molecular
numbering and include available atom, element, formal-charge, bond, and conformer context. Structured
diagnostic indices remain zero-based. Input parsing remains fail-fast.

Calculation provenance also includes `execution_metrics` with UTC start/end timestamps, phase durations for
parsing, applicability, computation, and non-JSON output writing, plus peak resident memory in MB.
Feature preparation is included in applicability timing; computation timing covers execution only.

## Install

Configure, build, and install to a chosen prefix:

```bash
cmake --preset gcc-release -DCMAKE_INSTALL_PREFIX="$PWD/_install"
cmake --build --preset gcc-release
cmake --install build/gcc-release --strip
_install/bin/chargefw calculate tests/fixtures/synthetic/sdf/water.sdf output
```

The installed directory can be moved after installation; bundled parameter JSON remains discoverable.
Installation includes the library, public headers, CLI, parameter data, and an exported CMake package:

```cmake
find_package(chargefw CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE chargefw::core)
```

The default installation includes nlohmann/json and Gemmi development packages and the private oneTBB
runtime, so adding the ChargeFW prefix to `CMAKE_PREFIX_PATH` is sufficient. A build configured with
`CHARGEFW_USE_SYSTEM_DEPENDENCIES=ON` instead requires those dependency prefixes to remain discoverable.

ChargeFW can also be added with `FetchContent`. Its CLI and tests default to off when it is a subproject;
installing the parent project installs ChargeFW's library, parameter JSON, and bundled dependencies into
the parent's prefix.

## Formatting

The optional pre-commit hook formats staged C++ files with the repository `.clang-format`:

```bash
pre-commit install
pre-commit run --all-files
```
