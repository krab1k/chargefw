# ChargeFW Project Guide

This document is the source of truth for ChargeFW's technical state, compatibility context, and
product roadmap. Read [AGENTS.md](AGENTS.md) for implementation rules, [TODO.md](TODO.md) for
actionable deliverables, and [README.md](README.md) for build and test commands.

## Purpose

ChargeFW is a C++23 framework for empirical partial atomic-charge calculation. It provides a
toolkit-neutral core model, built-in empirical methods, parameter-set loading and classification,
method applicability checks, and structured charge results.

It is being developed as a modern, library-first successor to **ChargeFW2**, the current
computational engine of [Atomic Charge Calculator III (ACC III)](https://acc.biodata.ceitec.cz).
ACC III is described in Raček *et al.*, *Atomic Charge Calculator III: a modern platform for
calculating partial atomic charges*, **Nucleic Acids Research** (2026),
DOI [10.1093/nar/gkag379](https://doi.org/10.1093/nar/gkag379).

ChargeFW is not yet the ACC III backend. `old/` is an archived ChargeFW2 copy used for behavior,
parameter, and compatibility research.

## Current architecture

```text
External molecular adapters / applications                 Not yet implemented
    (SDF, Mol2, PDB, mmCIF, SMILES, RDKit, Gemmi)
                              |
                              v
core::MoleculeCollection
    |- core::Molecule: atoms, bonds, conformers
    |- core::Atom: atomic number, formal charge, source name
    |- core::Bond: atom indices and order
    `- core::Conformer: coordinates sharing the molecule topology
                              |
                              v
features::PreparedMoleculeCollection
    |- PreparedMolecule: molecule + cached TopologyFeatures
    `- ConformerFeatures: on-demand geometry view per conformer
                              |
                 +------------+------------+
                 |                         |
                 v                         v
     parameters::ParameterSet       methods::MethodRegistry
                 |                         |
                 v                         v
      ParameterClassification     MethodRequirements/options
                 |                         |
                 +------------+------------+
                              v
                methods::find_applicable_methods()
                              |
                              v
                      methods::ApplicableMethod
                              |
                              v
                   methods::calculate_charges()
                              |
                              v
                      charges::ChargeSet
```

### Important public types

| Area | Types | Role |
|---|---|---|
| Core | `Atom`, `Bond`, `Conformer`, `Molecule`, `MoleculeCollection` | Input molecular graph and coordinates. |
| Features | `TopologyFeatures`, `ConformerFeatures`, `PreparedMolecule` | Cached/derived topology and geometry. |
| Parameters | `ParameterSet`, `ParameterClassification`, `ParameterView` | Parameter storage, matching, and method-facing lookup. |
| Methods | `Method`, `MethodRegistry`, `MethodRequirements`, `MethodOptions`, `ApplicableMethod` | Algorithm interface, capabilities, selection, and execution. |
| Charges | `AtomicCharges`, `ChargeAssignment`, `ChargeSet`, `ChargeCollection` | Atom-indexed calculated results and provenance. |

## Typical library workflow

```cpp
core::MoleculeCollection molecules{/* validated molecules */};
features::PreparedMoleculeCollection prepared{molecules};

const auto parameter_sets = parameters::load_default_parameter_sets();
const auto& registry = methods::method_registry();

std::vector<const methods::Method*> candidates;
for (const auto& method : registry.methods()) {
    candidates.push_back(method.get());
}

const auto applicability =
    methods::find_applicable_methods(prepared, candidates, parameter_sets);

// Application policy selects an explicitly reported candidate.
const auto& selected = applicability.applicable.front();
const charges::ChargeSet result = methods::calculate_charges(selected, prepared);
```

For normal application use, `calculation::calculate()` provides autodetection over the supplied
method and parameter candidates. It selects the applicable candidate with the highest method
priority, then the highest parameter-set priority. Ties are resolved deterministically by method ID
and parameter-set ID. Higher priorities therefore denote maintainer-curated automatic preference,
not a universal scientific quality ranking. The result retains applicability diagnostics when no
candidate can be calculated.

The current `chargefw` executable is a **water demonstration**. It builds two water conformers in
code, loads bundled parameter sets, autodetects the highest-priority applicable method and parameter
set, and prints charges. It is not yet a user-facing file/SMILES CLI.

### Calculation granularity

The smallest calculation unit is one molecule and zero or one selected conformer. A molecule
without a conformer remains valid for methods that do not require geometry; geometry-dependent
methods report it as inapplicable. Collection-wide and all-conformer execution are caller-level
iteration over these units. Future batch helpers must preserve molecule and conformer order and
report results for each unit independently.

## Implemented methods and parameters

The current registry contains 22 methods:

```text
abeem, charge2, delre, denr, dummy, eem, eqeq, eqeqc, formal, gdac,
kcm, mgc, mpeoe, peoe, qeq, sfkeem, smpqeq, sqe, sqeq0, sqeqp, tsef, veem
```

Bundled JSON parameter sets cover these parameterized methods and variants. `data/parameters/`
is installed under `share/chargefw/parameters`.

### Compatibility gap with ChargeFW2

The archived ChargeFW2 registry contains the same 22 methods. All nine archived SQE-family
parameter sets have been migrated. Numerical parity remains a release-blocking objective because
ACC III highlights SQE+qp.

## ChargeFW2 research summary

ChargeFW2 is a functional, application-oriented engine with CLI/Python bindings, custom
SDF/Mol2 reading, Gemmi PDB/mmCIF reading, output writers, OpenMP, and a nanoflann KD-tree. It
has useful production behavior, but combines concerns that are intentionally separated here:

| ChargeFW2 pattern | ChargeFW replacement |
|---|---|
| Atom stores graph data, coordinates, PDB residue data, Mol2 types, and classification state. | Core graph stays minimal; adapters and derived feature layers own source-specific metadata/caches. |
| Molecule owns all-pairs topology matrices and spatial KD-tree. | `TopologyFeatures` and `ConformerFeatures` own derived data. |
| Parameter classification mutates atoms and bonds. | Immutable `ParameterClassification` and `ParameterView`. |
| Method owns mutable parameters pointer and option values. | Stateless `Method::calculate(CalculationInput)`. |
| Suitability is mainly filtering/console behavior. | Structured applicability and prerequisite diagnostics. |
| PDB/mmCIF reader uses the first model only. | Core can represent multiple conformers; adapters should define model/altloc policy. |
| `full`, `cutoff`, and `cover` are coupled inside `EEMethod`. | Future execution policy should be explicit, reusable, and provenance-rich. |

The new implementation must preserve validated scientific behavior before making intentional
improvements. Compatibility fixtures should use identical topology, coordinates, formal charge,
options, and parameter data.

### SQE-family research

ChargeFW2 implements the SQE family with a signed bond-incidence matrix `T`. SQE solves
`(T A Tᵀ + diag(kappa)) p = T b`, then returns `q = Tᵀ p`, where `A` contains atom hardness on the
diagonal and the Gaussian-width Coulomb interaction `erf(d / sqrt(2 wi² + 2 wj²)) / d` off the
diagonal, and `b = -electronegativity`. SQE therefore conserves zero total charge within each
connected component. `sqeq0` adds formal charges as initial charges; `sqeqp` uses parameterized
`q0` after uniformly correcting its total to the molecular formal charge. All require coordinates,
atom parameters `electronegativity`, `hardness`, and `width`, plus bond parameter `kappa`; SQE
does not itself use formal charges, while SQE+q0 requires them. SQE+qp additionally requires atom
parameter `q0`, which it uniformly corrects to the molecule's formal-charge total. The archived
implementation does not explicitly diagnose singular or ill-conditioned transfer systems.

## Product direction

The release sequence is: establish scientific and parameter parity with ChargeFW2; expose a stable
request/result and provenance boundary; provide package, CLI, and JSON interfaces over that
boundary; then add adapters only when they serve an actual consumer. ACC III adoption remains
method-by-method and is gated by demonstrated parity.

The executable delivery checklist, validation criteria, and explicitly deferred work are maintained
in [TODO.md](TODO.md). In particular, cutoff/cover redesign, WebAssembly, MCP, ML/GNN adapters,
and custom parsers must not be pre-built ahead of validated need.

## Scientific and integration principles

- Preserve source atom indexing and map every output charge back to its source atom.
- Record method, parameter set, method options, conformer, approximation policy, and warnings in
  every serializable result.
- Keep exact/reference calculations distinct from cutoff/cover approximations.
- Make automatic selection explainable and overridable.
- Treat empirical charge results as method-dependent estimates; expose coverage and numerical
  limitations rather than implying universal validity.
