# ChargeFW

ChargeFW is a C++23, library-first framework for empirical partial atomic-charge calculation and a
modern successor to ChargeFW2, the engine used by Atomic Charge Calculator III (ACC III).

The repository includes a toolkit-neutral C++ library, 22 built-in methods, parameter loading and
classification, applicability/execution planning, full and parallel reduced calculation, molecular-file
adapters, and a focused CLI. It is not yet a production ACC III backend, Python package, or general
SMILES/chemistry-preparation tool.

## Documentation

| Document | Purpose |
| --- | --- |
| [PROJECT.md](PROJECT.md) | Implemented architecture, capabilities, limits, and product direction |
| [TODO.md](TODO.md) | Unfinished deliverables and acceptance criteria |
| [TESTING.md](TESTING.md) | Test-suite contracts, organization, and migration strategy |
| [AGENTS.md](AGENTS.md) | Repository boundaries and implementation rules |

## Requirements

- CMake 3.27 or newer
- Ninja
- GCC or Clang with C++23 support
- Internet access on first configure unless dependencies are already available to CMake

CMake searches for CLI11 2.7.2, nlohmann/json 3.12, Eigen 5.0, nanoflann 1.12, and Gemmi 0.7.4, then uses
`FetchContent` when needed. When `CHARGEFW_BUILD_TESTS` is enabled, it similarly searches for Snitch
1.3.2 and fetches it only if unavailable; production-only builds do not require Snitch.

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

The build-tree executable needs the parameter directory:

```bash
CHARGEFW_PARAMETER_DIR="$PWD/data/parameters" \
  build/gcc-debug/apps/chargefw/chargefw calculate \
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

Cutoff and explicit cover are available for ABEEM, EEM, EQeq, EQeq+C, QEq, SQE, SQE+q0, and SQE+qp.
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

`inspect` reports imported composition without loading parameters. `applicability` reports scientific
candidates and their full/cutoff/cover availability without calculation. `methods` and `parameters`
list installed capabilities.

### JSON result

Schema `1.0` contains one ordered result per imported molecule, source identity,
source-ordered charge assignments, totals, and diagnostics. `calculation_provenance` records requested
conformer selection, selection/classification/resource/structural policies, and the effective method,
parameter set, method options, execution mode, radius, correction, and warnings. Requested method
options contain explicit overrides; effective options contain the complete selected-method values.
JSON charge values are rounded to at most four decimal places; internal calculations retain native
precision.

Calculation provenance also includes `execution_metrics` with UTC start/end timestamps, phase durations for
parsing, applicability, computation, and non-JSON output writing, plus peak resident memory in MB.
Feature preparation is included in applicability timing; computation timing covers execution only.

## Install

After configuring and building, install to a chosen prefix:

```bash
cmake --install build/gcc-release --prefix "$PWD/_install" --strip
env -u CHARGEFW_PARAMETER_DIR \
  _install/bin/chargefw calculate tests/fixtures/synthetic/sdf/water.sdf output
```

Installation includes the library, public headers, CLI, and parameter data. Exported CMake package
targets and binary distribution packages are not implemented yet.

## Formatting

The optional pre-commit hook formats staged C++ files with the repository `.clang-format`:

```bash
pre-commit install
pre-commit run --all-files
```
