# ChargeFW

ChargeFW is a C++23, library-first framework for empirical partial atomic-charge calculation and a
modern successor to ChargeFW2, the engine used by Atomic Charge Calculator III (ACC III).

The repository includes a toolkit-neutral C++ library, 22 built-in methods, parameter loading and
classification, applicability/execution planning, full and serial cutoff calculation, molecular-file
adapters, and a focused CLI. It is not yet a production ACC III backend, Python package, or general
SMILES/chemistry-preparation tool.

## Documentation

| Document | Purpose |
| --- | --- |
| [PROJECT.md](PROJECT.md) | Implemented architecture, capabilities, limits, and product direction |
| [TODO.md](TODO.md) | Unfinished deliverables and acceptance criteria |
| [AGENTS.md](AGENTS.md) | Repository boundaries and implementation rules |
| [CUTOFF_IMPLEMENTATION_PLAN.md](CUTOFF_IMPLEMENTATION_PLAN.md) | Reduced-execution design record and remaining validation |

## Requirements

- CMake 3.27 or newer
- Ninja
- GCC or Clang with C++23 support
- Internet access on first configure unless dependencies are already available to CMake

CMake searches for CLI11 2.6, nlohmann/json 3.12, Eigen 5.0, nanoflann 1.12, and Gemmi 0.7.4, then uses
`FetchContent` when needed.

## Build and test

Run from the repository root:

```bash
cmake --preset gcc-debug
cmake --build --preset gcc-debug
ctest --preset gcc-debug
```

Run one focused test after building:

```bash
ctest --test-dir build/gcc-debug -R test_cutoff_execution --output-on-failure
```

Sanitizer configuration:

```bash
cmake --preset clang-asan
cmake --build --preset clang-asan
ctest --preset clang-asan
```

Other tracked configure presets are `gcc-release`, `clang-debug`, and `clang-tidy`. Useful options
include `CHARGEFW_BUILD_TESTS`, `CHARGEFW_BUILD_CLI`, `CHARGEFW_ENABLE_CCACHE`,
`CHARGEFW_ENABLE_SANITIZERS`, and `CHARGEFW_ENABLE_CLANG_TIDY`.

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
- `--full-atom-threshold COUNT|unlimited` changes the automatic full-execution safeguard (default
  20,000 atoms).

Automatic selection uses full execution when it is not discouraged. For an expensive candidate above
the threshold, it uses supported cutoff at 12 Å before considering another candidate policy. Explicit
full overrides the threshold and records a warning. Cover is accepted by the policy API but no method
currently supports it.

Cutoff is currently available for ABEEM, EEM, EQeq, EQeq+C, QEq, SQE, SQE+q0, and SQE+qp.

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

Schema `1.0` contains one ordered result per imported molecule, source identity and atom mapping,
source-ordered charge assignments, totals, and diagnostics. `calculation_provenance` records requested
conformer selection, selection/classification/resource/structural policies, and the effective method,
parameter set, execution mode, radius, correction, and warnings. JSON charge values are rounded to at
most four decimal places; internal calculations retain native precision.

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
