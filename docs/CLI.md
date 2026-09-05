# Command-line interface

The `chargefw` executable imports molecular records, assesses available methods and parameter sets, and
writes source-mapped partial charges.

## Installation

ChargeFW requires CMake 3.28 or newer, Ninja, and a GCC or Clang toolchain with C++23 support. The default
build downloads pinned dependencies when they are not already available.

```bash
cmake --preset gcc-release -DCMAKE_INSTALL_PREFIX="$PWD/_install"
cmake --build --preset gcc-release
cmake --install build/gcc-release --strip
```

Use the installed executable. Bundled parameter data is discovered relative to the installed library:

```bash
_install/bin/chargefw --help
```

## Commands

```text
chargefw calculate [options] INPUT OUTPUT_DIRECTORY
chargefw inspect [input-options] INPUT
chargefw applicability [options] INPUT
chargefw methods
chargefw parameters [METHOD]
```

Run `chargefw COMMAND --help` for the option syntax accepted by the installed version.

## Calculate charges

```bash
chargefw calculate molecule.sdf output
```

If no method or parameter set is specified, ChargeFW assesses the bundled catalog and chooses the first
plan in deterministic priority order. The output directory is created if necessary.

### Input options

Input format is selected from the file extension.

| Extensions | Reader |
| --- | --- |
| `.mol`, `.sdf` | Native MOL/SDF reader |
| `.mol2` | Native Tripos MOL2 reader |
| `.json` | ChargeFW molecule JSON 1.0 reader |
| `.pdb` | Gemmi-backed PDB reader |
| `.cif`, `.mmcif` | Gemmi-backed mmCIF reader |

The [molecular format reference](FORMATS.md) describes the supported subsets, imported molecular data,
record mapping, warnings, and errors.

`--conformers first|all` selects the first conformer/model or all conformers/models. The default is
`all`.

PDB and mmCIF input additionally accepts:

- `--structural-selection all|polymers-and-ligands|polymers` (default `all`);
- `--structural-bonds none|explicit|templates|hybrid` (CLI default `hybrid`).

Structural options are rejected for other input formats. Alternate locations prefer blank, then `A`,
then the first available location.

The CLI reads the imported collection before calculation and stops at the first malformed record.

### Method and execution options

| Option | Meaning |
| --- | --- |
| `--method ID` | Restrict assessment to one method |
| `--parameter-set ID` | Restrict assessment to one bundled parameter set |
| `--permissive-types` | Allow [permissive parameter classification](PARAMETERS.md#strict-and-permissive-matching) |
| `--method-option METHOD.OPTION=VALUE` | Override a method option; repeatable |
| `--execution auto|full|cutoff|cover` | Select execution policy; default `auto` |
| `--radius ANGSTROM` | Radius for cutoff/cover; explicit reduced modes require at least 8 Å |
| `--charge-correction uniform|none` | Correction for explicit reduced execution |
| `--cutoff-atom-threshold COUNT|unlimited` | Automatic full-to-cutoff threshold; default 20,000 |
| `--cover-atom-threshold COUNT|unlimited` | Automatic cutoff-to-cover threshold; default 80,000 |
| `--threads COUNT` | Maximum calculation threads; `0` delegates to oneTBB |
| `--progress` | Render calculation progress on standard error |

Use `chargefw methods` to list option IDs, defaults, allowed choices, and numeric bounds. Method options
are method-scoped even when a method is selected explicitly:

```bash
chargefw calculate --method peoe --method-option peoe.iters=8 molecule.sdf output
```

Explicit method, parameter-set, or execution choices do not fall back to alternatives if they are
inapplicable. `applicability` reports the reasons without running a calculation.

### Automatic execution

Automatic planning prefers full execution. For methods with cubic-time or quadratic-memory behavior,
collections containing a molecule above the cutoff threshold use supported cutoff execution. Above the
cover threshold, supported cover execution is preferred. Automatic reduced execution uses a 12 Å radius
and uniform charge correction.

Explicit `full`, `cutoff`, or `cover` filters plans to that mode. Explicit full execution can exceed the
resource threshold and reports a warning rather than silently changing mode.

Cutoff and cover are available for ABEEM, EEM, EQeq, EQeq+C, QEq, SFKEEM, SQE, SQE+q0, and SQE+qp.
They are explicit approximations and do not have a general accuracy guarantee.

## Output

The output basename is derived from the input filename:

```text
OUTPUT_DIRECTORY/<input-stem>.chargefw
```

On a successful nonstructural calculation, ChargeFW writes:

```text
<basename>.json
<basename>.sdf
<basename>.mol2
<basename>.cif
```

PDB and mmCIF input writes JSON and mmCIF only. SDF and MOL2 output preserves source content when the
input has the same format. Generated SDF and MOL2 use the first retained conformer; generated mmCIF writes
all retained conformers. The [molecular format reference](FORMATS.md#charge-output) describes
preservation, generated structures, charge fields, precision, and mapping checks.

Generated SDF and MOL2 use source conformer zero and therefore require source coordinates even when the
selected method is geometry-independent. With `--conformers first`, preservation-oriented structural
output may retain uncalculated source models, but charges are written only for the selected conformer.

Once import and request construction succeed, normal calculation outcomes write `<basename>.json`.
Molecular charge files are written only on success. Import, request-construction, filesystem, and output
compatibility failures are reported on standard error and can occur before a result document is written.

During calculation, press `Ctrl+C` once to request cooperative cancellation. ChargeFW stops at its next
cancellation check point, writes `<basename>.json` with status `cancelled` and no charge assignments, and
exits with status 5. This does not interrupt import, request construction, or output writing.

### JSON result

ChargeFW writes result JSON schema `1.0`; this is distinct from molecule input JSON schema `1.0`. See the
[result JSON format](FORMATS.md#chargefw-result-json-10) for its records, diagnostics, assignments,
provenance, metrics, and numeric precision.

## Inspection and discovery

Inspect imported records without loading parameter sets:

```bash
chargefw inspect molecule.sdf
```

Report runnable plans, rejected alternatives, and the selected plan without calculating:

```bash
chargefw applicability --method eem molecule.sdf
```

List methods and their option schemas, or list bundled parameter sets:

```bash
chargefw methods
chargefw parameters
chargefw parameters eem
```

## Exit statuses

| Status | Meaning |
| --- | --- |
| `0` | Calculation or reporting command completed successfully |
| `1` | Unexpected internal failure |
| `2` | Invalid input or request |
| `3` | No executable plan |
| `4` | Numerical calculation failure |
| `5` | Calculation cancelled |
