# Command-line interface

The `chargefw` executable imports molecular records, assesses available methods and parameter sets, and
writes source-mapped partial charges.

## Installation

Follow the [native source installation](NATIVE.md#build-and-link) with
`CHARGEFW_BUILD_CLI=ON`. The installed executable discovers bundled parameter data relative to the
installed library.

```bash
_install/bin/chargefw --help
```

The root `Dockerfile` builds the same CLI as a portable Ubuntu image without Python, tests, or
host-specific processor instructions:

```bash
docker build --tag chargefw:local .
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/work" -w /work chargefw:local \
    calculate molecule.sdf output
```

## Commands

```text
chargefw calculate [options] INPUT OUTPUT_DIRECTORY
chargefw inspect [input-options] INPUT
chargefw applicability [options] INPUT
chargefw methods [METHOD]
chargefw parameters [--method METHOD] [PARAMETER_SET]
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

Structural options are rejected for other input formats. Selection, bonding, and alternate-location
semantics are defined in the [PDB and mmCIF format reference](FORMATS.md#pdb-and-mmcif-input).

The CLI reads the imported collection before calculation and stops at the first malformed record.

### Method and execution options

| Option | Meaning |
| --- | --- |
| `--method ID` | Restrict assessment to one method |
| `--parameter-set ID` | Parameter set for `--method` (required) |
| `--permissive-types` | Allow [permissive parameter classification](PARAMETERS.md#strict-and-permissive-matching) |
| `--method-option METHOD.OPTION=VALUE` | Override a method option; repeatable |
| `--execution auto|full|cutoff|cover` | Select execution policy; default `auto` |
| `--radius ANGSTROM` | Radius for cutoff/cover; explicit reduced modes require at least 8 Å |
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

Automatic and explicit execution share the policy, supported-method, conservation, and approximation
semantics documented in [Assessment and execution](PROJECT.md#assessment-and-execution). Explicit choices
do not silently fall back to another mode.

## Output

The output basename is derived from the input filename:

```text
OUTPUT_DIRECTORY/<input-stem>.chargefw
```

Every normal calculation outcome writes result JSON:

```text
<basename>.json
```

On success, request fresh molecular representations explicitly:

```bash
chargefw calculate --output-mol2 molecule.sdf output
chargefw calculate --output-mmcif structure.cif output
chargefw calculate --output-mol2 --output-mmcif molecule.sdf output
```

The requested files use the same basename:

```text
<basename>.mol2
<basename>.cif
```

MOL2 and mmCIF are fresh generated representations for every supported input format and are not written
unless requested. The [molecular format reference](FORMATS.md#charge-output) describes their intended use,
coordinate requirements, generated structures, charge fields, precision, and mapping checks.

Once import and request construction succeed, normal calculation outcomes write `<basename>.json`. On
success, ChargeFW writes JSON before attempting requested molecular exports. A molecular export failure
reports `Export error`, exits with status 6, and leaves JSON and any earlier completed export intact; the
JSON calculation status remains `success`. Files are written directly to their destinations, so a write
failure can leave that destination incomplete. Import, request-construction, primary JSON output, and other
filesystem failures are reported on standard error and can occur before a result document is published.

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

List method IDs and full names, inspect one method, or list bundled parameter sets:

```bash
chargefw methods
chargefw methods eem
chargefw parameters
chargefw parameters --method eem
chargefw parameters QEq_original
```

`chargefw methods METHOD` reports the selected method's identity, publication, human-readable notes,
priority, coordinate requirement, time and memory complexity, cutoff and cover support, and option schema.
Complexity uses Big-O notation with `n` for atoms and `m` for bonds, as defined in the
[project design](PROJECT.md#methods-and-parameters).

`chargefw parameters` lists parameter-set IDs and names. Use `--method METHOD` to filter the summary by
the parameter set's declared method. Pass a parameter-set ID instead to report its method, name,
publication, notes, and priority; the positional ID and `--method` cannot be combined.

## Exit statuses

| Status | Meaning |
| --- | --- |
| `0` | Calculation or reporting command completed successfully |
| `1` | Unexpected internal failure |
| `2` | Invalid input or request, or primary output publication failure |
| `3` | No executable plan |
| `4` | Numerical calculation failure |
| `5` | Calculation cancelled |
| `6` | A requested molecular export failed after JSON publication |
