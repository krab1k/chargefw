# Project design

## Purpose

ChargeFW calculates empirical partial atomic charges from an explicit molecular graph, formal charges,
and, when required by a method, Cartesian coordinates. The same native engine is used by the C++ API,
the CLI, and the Python package.

ChargeFW is library-first. File readers and toolkit adapters translate external representations into an
owned, toolkit-neutral molecule model; method applicability, parameter matching, execution selection,
and calculation remain in the native library.

ChargeFW implements finite, non-periodic molecular calculations. Periodic and Ewald variants are not
implemented. ChargeFW is not yet the production backend for Atomic Charge Calculator III.

## Design principles

- **Explicit molecular input.** ChargeFW does not silently change topology, protonation, formal charge,
  hydrogen count, coordinates, method options, parameter sets, or execution mode.
- **Scientific applicability is separate from execution policy.** Missing coordinates, topology,
  element properties, or parameter coverage are hard applicability failures. Resource thresholds only
  guide automatic execution.
- **Selection is deterministic.** Unless constrained by the caller, methods and parameter sets are
  ordered by priority and then by stable ID.
- **Full and reduced execution are distinct.** Full execution follows the implemented finite-molecule
  method. Cutoff and cover are explicit approximations and are recorded in result provenance.
- **Mappings survive every boundary.** Molecule, atom, and conformer order are preserved in calculation
  results and adapters retain available source identities.
- **Algorithms are reusable and concurrent.** Prepared assessments and concrete plans own the state they
  need and can be executed repeatedly without global mutable calculation state.

## Architecture

```text
files, C++ values, or Python arrays
                 |
                 v
       toolkit-neutral molecules
       atoms + bonds + conformers
                 |
                 v
       prepared topology/geometry
                 |
        +--------+--------+
        |                 |
        v                 v
 bundled parameters   method registry
        |                 |
        +--------+--------+
                 v
       applicability assessment
                 |
                 v
       concrete execution plans
                 |
                 v
          full/cutoff/cover
                 |
                 v
      source-ordered charge sets
```

The public native namespaces mirror these responsibilities:

| Namespace | Responsibility |
| --- | --- |
| `chargefw::core` | Owned atoms, bonds, molecules, conformers, and periodic-table data |
| `chargefw::features` | Prepared topology, geometry, and spatial fragments |
| `chargefw::parameters` | Parameter data, loading, and atom/bond classification |
| `chargefw::methods` | Method metadata, options, requirements, applicability, and algorithms |
| `chargefw::calculation` | Assessment, deterministic planning, execution policy, progress, and calculation |
| `chargefw::charges` | Source-indexed charge assignments and charge sets |
| `chargefw::adapters` | Native and Gemmi-backed molecular input and output |

Applications normally use the owned assessment facade rather than the lower-level prepared-feature and
classification interfaces.

## Methods and parameters

The built-in registry contains 22 method IDs:

```text
abeem, charge2, delre, denr, dummy, eem, eqeq, eqeqc, formal, gdac,
kcm, mgc, mpeoe, peoe, qeq, sfkeem, smpqeq, sqe, sqeq0, sqeqp, tsef, veem
```

Bundled JSON parameter sets are installed with the library. Parameter-dependent methods are assessed
against those sets unless a native caller supplies another catalog. The CLI and Python package expose
the bundled catalog only. The [parameter-set JSON reference](PARAMETERS.md) defines custom native
parameter data and its classifier semantics.

Methods declare their coordinate, topology, element-property, formal-charge, and parameter
requirements. Assessment reports both runnable plans and structured reasons for rejected candidates.
Method options are validated against each method's schema before execution.

## Assessment and execution

Assessment prepares a molecule collection once, performs parameter classification, evaluates scientific
requirements, and expands applicable candidates into concrete plans. A plan contains the selected
method, parameter set, validated options, execution policy, and any policy warnings.

The effective execution mode is always one of:

- `full`: calculate the complete molecular target;
- `cutoff`: solve one radius fragment per source atom and retain its mapped center charge;
- `cover`: solve radius fragments around source-order pivots and retain charges in covered interiors.

Automatic planning prefers full execution. For methods classified as expensive, the default resource
policy moves to cutoff above 20,000 atoms and to cover above 80,000 atoms when those modes are supported.
Automatic reduced execution uses a 12 Å radius. Explicit reduced execution requires a radius of at least
8 Å. Explicit full execution may override a resource threshold and returns a warning rather than changing
the requested mode.

Cutoff and cover are implemented for:

```text
abeem, eem, eqeq, eqeqc, qeq, sfkeem, sqe, sqeq0, sqeqp
```

Reduced execution defaults to uniform correction of the final molecular charge. A caller can explicitly
request no correction where the execution policy permits it. Reduced calculations are approximations;
the project does not claim a general accuracy envelope for them.

## Results and provenance

Geometry-dependent methods produce one assignment for each molecule conformer. Geometry-independent
methods produce one assignment per molecule. Every assignment identifies its source molecule and, when
applicable, conformer; atom values remain in source order.

Application-facing results distinguish success, invalid input or request, no executable plan, numerical
failure, and cancellation. Successful results retain the effective method, parameter set, validated
options, execution mode, radius, correction policy, and execution warnings. The CLI additionally records
requested policy, diagnostics, timings, timestamps, and peak resident memory in its JSON result.

## Molecular data scope

The native adapters support MOL/SDF, MOL2, ChargeFW molecule JSON 1.0, PDB, and mmCIF input and write
ChargeFW result JSON 1.0, SDF, MOL2, and mmCIF charge data. The language-independent
[molecular format reference](FORMATS.md) defines exactly what each reader imports and how each writer
preserves or generates molecular data.

ChargeFW is not a general chemistry-preparation toolkit. It does not provide SMILES parsing, arbitrary
bond perception, sanitization, protonation, hydrogen addition or removal, conformer generation, or
geometry optimization. Structural template bonding covers common amino acids, nucleotides, water, and
basic polymer links; it is not a complete Chemical Component Dictionary provider.

Current release and integration gaps are listed in the root [TODO](../TODO.md).
