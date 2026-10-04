# Fixed-Charge Electrostatic Embedding: Research and Implementation Handoff

## Status and Purpose

Research date: **2026-10-03**.

This is a research/design handoff requested by the user, not documentation of an implemented feature.
No embedding capability, public API, parameter set, or charge solver was changed during this research.
The proposed names below are illustrative, not committed API names. This document records the evidence,
scientific rationale, design direction, unresolved decisions, and checks needed to resume independently.

Repository revision inspected: `05bbd4139d80b2fd70ab70728047ced015aea09d`.

Follow [AGENTS.md](AGENTS.md) before implementation. The user-requested [PLAN.md](PLAN.md) now owns this
feature's staged implementation checklist and model responsibilities: Astra handles architecture,
critical scientific/API decisions, and review; Luna handles coding within the agreed design.
Keep broader unfinished work in [TODO.md](TODO.md); once behavior exists, document it in the appropriate
owning documents under `docs/`. Keep this research history separate from those implemented-behavior contracts. Build and
validation commands belong to [DEVELOPMENT.md](DEVELOPMENT.md).

## Executive Summary

The original problem was that all three bundled SQE+qp parameter sets lack Mg. This prevents their use
on otherwise suitable protein structures containing Mg, even when the desired output is primarily the
protein's charges.

The recommended direction is **explicit, capability-checked fixed-charge electrostatic embedding**:

- Calculate charges on an active subsystem, prioritizing EEM, the SQE family, and QEq.
- Represent selected ions, and eventually selected molecular fragments, by prescribed atomic charges.
- Include their electrostatic potential in the active solve, so protein charges respond to them.
- Keep environmental charge values fixed and preserve all original output mappings.
- Require parameters for active atoms/bonds only, without weakening the parameter matcher.

The user explicitly rejected a purely decoupled approximation that removes ions, calculates the protein,
and appends ion charges afterward. **Direct protein response to the environment is essential.**

This has a sound methodological precedent in electrostatically embedded QM/MM and a related precedent
in electronegativity-based fluctuating-charge MD. It is not a validated ChargeFW model yet. Existing
parameters fitted to molecular charges are not automatically validated for response to ionic fields.

The scientific opportunity is more meaningful than an archive-wide percentage suggests:

- Important represented classes include calcium sensors, adhesion proteins, channels, phosphatases,
  and magnesium-dependent enzymes.
- Live PDBe research found 6,380 conservative protein-plus-common-ion candidates, ignoring water.
- Local CCD/native assessment confirmed coverage of several noncanonical residues and of explicit
  ATP4-, ADP3-, and GTP4- graphs in existing broader parameter sets.
- A more inclusive, still restricted PDBe query allowing seven already parameter-matched ligands found
  11,465 ion-containing composition candidates, including 1,951 nucleotide-containing entries.

These are composition and parameter-coverage findings, **not successful whole-protein calculations or
charge-accuracy measurements**.

Follow-up discussion on 2026-10-03 narrowed the relevant methods to EEM, SQE/SQE+q0/SQE+qp, QEq,
ABEEM, EQeq, and possibly EQeq+C. The central requirement is to avoid needing fitted ion response
parameters while retaining the active molecule's response to prescribed ion charges. Keep the first
implementation small: molecular frozen-source support is optional, not a prerequisite. Sections 3.3,
4.6, and 11 distinguish the resulting recommendations from choices still awaiting a concrete contract.

### Navigation

- [Motivation](#1-why-not-just-add-mg-parameters)
- [Scientific precedents and limitations](#2-scientific-meaning-and-precedents)
- [Ions and molecular fragments](#3-fixed-ions-and-frozen-molecular-fragments)
- [Equations](#4-equations-and-implementation-specific-conventions)
- [Proposed architecture](#5-proposed-application-architecture)
- [Implementation code map](#6-code-map-for-resumption)
- [PDBe coverage and reproduction](#7-live-pdbe-coverage-research)
- [Biological significance](#8-biological-significance-of-the-candidate-pool)
- [CCD/native experiments and reproduction](#9-local-ccd-and-native-parameter-coverage-experiments)
- [Implementation and validation sequence](#10-validation-and-implementation-sequence)
- [Open decisions](#11-open-decisions-and-recommendation)
- [Bibliography](#references)

## 1. Why Not Just Add Mg Parameters?

The source SQE+qp publication [1] fitted organic molecules and small peptides, not comprehensive
metalloprotein chemistry. Its peptide dataset contains 60 di- and tripeptides and no Mg.

| Bundled SQE+qp set | Atom classifier | Element coverage | Reference charges |
| --- | --- | --- | --- |
| `SQEqp_Schindler2021_PUB_pept` | `bonded` | C, H, N, O, S | B3LYP/6-31G*/NPA |
| `SQEqp_Schindler2021_DTP_small` | `hbo` | C, H, N, O, S | B3LYP/6-311G/NPA |
| `SQEqp_Schindler2021_CCD_gen` | `hbo` | Br, C, Cl, F, H, N, O, P, S | B3LYP/6-311G/NPA |

SQE+qp requires atom `electronegativity`, `hardness`, `width`, and `q0`, plus bond `kappa`.
Here `q0` is a fitted initial charge, not the supplied formal charge. Atomic ionization energies and
electron affinities can motivate initial guesses for some terms, but do not determine a compatible
width, fitted seed, and metal-ligand bond hardness. Neutral Mg properties are also not a sufficient
description of coordinated Mg(II).

The fitted quantities should not all be interpreted as independently measurable atomic properties:
some existing widths and hardnesses are negative, widths enter squared, and a common offset in all
SQE+qp seed parameters disappears on normalization. In particular, some DTP seed values near -23 are
not literal atomic charges.

Adding Mg-O graph bonds changes the model, enables charge transfer, and changes `bonded` ligand atom
types and adjacent bond keys. A reliable fully variable-charge Mg extension would need representative
QM coordination data and fitting, not just one extra atom row. Such a study remains a possible future
direction, but is not required to implement the proposed fixed-source approximation.

Existing alternatives do not remove this distinction:

- `eqeq` can use built-in Mg elemental data, but its neutral-reference-center variant and broad coverage
  do not establish accuracy for Mg2+ in proteins.
- `EQeqC_original_MOF` contains Mg, but is a different model/REPEAT fitting target and lacks S and P.
  Its numbers should not be transplanted into SQE+qp.
- `formal` returns supplied formal charges, but does not calculate the requested protein response.

## 2. Scientific Meaning and Precedents

### 2.1 Active and frozen subsystems

The active subsystem has variable charges. The environment supplies fixed charge values and coordinates.
This includes active-region charge redistribution, but excludes mutual electronic polarization and
charge transfer across the active/frozen boundary during the solve.

"Fixed" describes charge values. It does not require that coordinates be identical across conformers
or an MD trajectory. "One-way" describes electronic response, not forces acting only in one direction.
ChargeFW remains a charge calculator, not a validated MD force field or force/trajectory engine.

| Treatment | Active response | Environmental response |
| --- | --- | --- |
| Remove environment, calculate, append charges | No response to omitted field | None |
| Conventional additive fixed-charge MD | Protein charge values remain fixed | Charge values remain fixed |
| Electrostatic QM/MM embedding | QM density responds to MM charges | MM charges remain fixed |
| Proposed EEM/SQE embedding | Active atomic charges respond to prescribed field | Source charges remain fixed |
| Mutually polarizable embedding | Active subsystem responds | Environment responds too |

Electrostatic QM/MM is the closest conceptual analogue [2-5]. Fluctuating-charge MD based on
electronegativity equalization has existed since at least Rick, Stuart, and Berne (1994) [6]. Fixed
ion and molecular partial-charge models are also established in additive MD [7-9]. Coarse-grained
charged sites are a looser analogy, not the main scientific justification for this feature.

These references establish a methodological foundation, not the quantitative accuracy of ChargeFW's
selected parameterization. In particular, the QEq model in [3] describes static charge distribution;
that work uses a separate induced-dipole model for response. Do not cite it as validating arbitrary
EEM/QEq response parameters.

### 2.2 Limitations to validate

- Fixed source charges omit source polarization and dynamic inter-subsystem charge transfer.
- A fixed Mg2+ charge is a model choice; it need not equal a coordinated Mg NPA partial charge.
- Short-range field errors, charge penetration, and missing electronic effects matter near contacts.
  Reference [4] specifically studies water and Mg2+-water embedding examples.
- Divalent-ion MD modeling itself has known limitations. Improved models such as 12-6-4 account for
  effects missing from simple 12-6 models [8]. Those corrections cannot simply be transplanted into
  a responsive EEM/SQE calculation without considering the model and possible double counting.
- Successful static-charge fitting does not establish correct susceptibility to an external field.
- EEM's existing long-range charge-transfer behavior is not fixed by embedding.
- A geometry-dependent charge response is not a prediction of binding affinity, catalysis, allostery,
  redox chemistry, or proton transfer.
- Adding explicit ions is not equivalent to adding solvent screening or an ionic-strength model.
  Do not introduce an unexplained dielectric factor.

Recommended description, once the feature exists:

> We use a fixed-charge electrostatic-embedding approximation, conceptually analogous to electrostatic
> embedding in QM/MM. Prescribed ion charges or molecular partial charges generate an external potential
> to which the active subsystem responds through an EEM or SQE charge-equilibration model. The approach
> includes active-region charge redistribution but excludes mutual polarization and charge transfer
> between the active and frozen regions.

Avoid "the same approximation as standard fixed-charge MD" and "validated because it is standard in
MD". MD ion models also include short-range interactions and solvent-specific calibration.

## 3. Fixed Ions and Frozen Molecular Fragments

### 3.1 Monatomic ions

A prescribed point Mg2+ source needs its position and charge. For the proposed EEM and SQE couplings,
no Mg electronegativity, self hardness, fitted seed, or Mg-O bond parameter is needed.

A parameter-free point-source model contains no ion-specific radius or polarizability: Na+ and K+ at
identical positions generate the same prescribed field, as do Mg2+ and Ca2+. Differences can arise from
their observed geometry, but this is not a model for predicting ion selectivity or hydration structure.

Initial scientific emphasis should be Na+, K+, Mg2+, Ca2+, and Cl-. Transition-metal, highly charged,
and strongly covalent sites can be supplied as fields mechanically but require stronger caveats and
separate validation. Technical capability is not a claim of universal metal-site accuracy.

### 3.2 Molecular environments, including sulfate

For SO4(2-), a total charge of -2 does not determine the near-field electrostatic potential. Even with
four equivalent oxygens, `q_S + 4*q_O = -2` leaves the charge distribution undetermined.

The sensible extension is a prescribed **per-atom partial-charge vector**, not automatically a -2 point
at sulfur or the molecular centroid. A single monopole can be a far-field approximation, not the default
for a bound sulfate. Constraining only a ligand's total charge while allowing its atoms to respond is a
different model and still requires response parameters.

Possible source-charge provenance:

- An explicitly chosen, curated fixed-charge model for the component.
- Suitable QM-derived partial charges, with charge scheme, geometry/environment, and protonation stated.
- A separate, explicit calculation with an existing supported model, followed by freezing its result.

CCD formal charges establish chemical identity and net charge; they are not generally a physical
partial-charge distribution for a multiatom fragment. Using a localized sulfate Lewis formal-charge
pattern directly can introduce arbitrary oxygen asymmetry. SO4(2-) and HSO4(-) are different sources.

If protein charges are the observable of interest, the environmental charge distribution still matters:
errors in its field propagate into the protein. Frozen fragments are an explicit approximation, not a
way to make ligand chemistry irrelevant. Reference [9] provides related classical sulfate/phosphate
modeling precedent, not ready-made validation for ChargeFW embedding.

Start with whole components having no covalent bonds crossing the active/frozen boundary. Cutting a
covalent bond introduces a separate boundary/capping problem. Bonds internal to a wholly frozen
fragment need not become active charge-transfer variables.

### 3.3 Follow-up: fragment charges and minimal scope

The unresolved fragment question is primarily **where to obtain per-atom source charges**, not how to
sum their electrostatic contributions. A monatomic ion's prescribed net charge determines its single
atomic charge. A molecular component's net charge does not determine its atomic distribution.

PDB/mmCIF component-instance selection can use the same mapping machinery for MG and SO4. However, a
CCD component instance is not necessarily a disconnected component of the imported bond graph. A ligand,
residue, or ion may have graph bonds to the active region. For a simple initial boundary, reject any
active/fixed crossing bond rather than silently deleting it; internal bonds of a wholly frozen component
do not enter the active solve. This also requires care when coordination bonds are imported.

Recommended scope, without making optional fragment features a release requirement:

- Prioritize explicitly prescribed monatomic ions. An ion-only first implementation is acceptable.
- Keep parameter-covered molecules such as sulfate active unless the user explicitly wants them frozen.
  Section 9.4 records sulfate coverage in the CCD SQE+qp and broader EEM sets, not in every set and not
  proof of accuracy. One selected set must cover the complete active graph; do not splice fitted rows.
- Represent fixed sources simply as mapped original atom indices and prescribed finite real charges,
  with coordinates obtained from their conformers. Avoid element-specific solver branches or an
  integer-only ion representation. A frozen molecule can later use the same sum of atomic sources.
- If explicit molecular sources are exposed, require a complete mapped charge vector and validate its
  declared total and graph boundary. This adds preparation/interface checks, not a new charge solver.
- Defer automatic fragment charge generation, a general functional-group rule engine, molecular preset
  libraries, covalent cutting/capping, and independent external charge clouds. Start with sources that
  belong to the original target. For any external source field, the prepared active charge remains the
  formal-charge sum of active atoms and prescribed source charges contribute separately to the modeled
  total, as in section 4.1.

Formal-charge or chemistry-based rules are possible **explicit approximations**, not forbidden inputs.
For example, `q_S = 0` and four `q_O = -0.5` give a symmetric sulfate source with total -2, avoiding an
arbitrary choice of two localized -1 oxygens. This is illustrative, not an adopted sulfate preset:
positive sulfur with more negative oxygens can have the same total but a different near-field potential.
Any supplied rule needs chemical-state matching, exact charge accounting, and named provenance. Chemical
intuition and symmetry alone do not establish field accuracy. No fragment charge rule was selected.

The user's recollection of Open Babel charge seeds has a concrete source: the `InitialPartialCharges()`
helper in [Open Babel molchrg.cpp](https://github.com/openbabel/openbabel/blob/master/src/molchrg.cpp),
inspected on 2026-10-03, assigns -0.5 to atoms matching `IsCarboxylOxygen()`, -0.666 to selected phosphate
oxygens, -0.5 to sulfate oxygens, and formal charges otherwise. Its call inside `AssignPartialCharges()`
is commented out in the inspected source; that routine initializes from existing atomic partial charges.
This is evidence for rule-based seeding, not verification of the active default preparation path or a
validated final frozen-charge model. The linked branch is mutable; recheck before relying on it.

A separately chosen charge calculation followed by freezing is another source option. PEOE could serve
as such a generator without being an active embedding method. A seed, a final calculated partial charge,
and a validated embedding-source charge are different things. Neither this workflow nor new presets
need to be built into the initial embedding feature.

## 4. Equations and Implementation-Specific Conventions

### 4.1 Charge budget and fixed-value substitution

Let A denote active atoms, F fixed atoms belonging to the original target, and f their prescribed charges.
The original proposal used the original-target convention below; that policy is retained here as research
history and was superseded for the current implementation by the dated follow-up after this section.

```text
Q_A = Q_total - sum(f)
phi_i = sum_F G_iF(r_iF) * f_F
```

For a quadratic atom-charge model with energy `chi^T q + 0.5*q^T H*q`, fixing the environmental charges
gives:

```text
[ H_AA   1 ] [ q_A    ] = [ -chi_A - H_AF*f ]
[ 1^T    0 ] [ lambda ]   [ Q_A             ]
```

The external potential is `phi = H_AF*f`. It is subtracted from the current charge RHS, or equivalently
added to electronegativity. This is fixed-value substitution, not a Schur complement that lets sources
respond. Source self energies and source-source interactions are constant with respect to active charges.
The cross interaction still has to be defined; eliminating a self parameter does not necessarily remove
that parameter from a method's pair kernel.

Solving the whole system normally and overwriting fixed atoms afterward is not equivalent.

#### Active-charge policy follow-up (2026-10-04)

The current ChargeFW model instead takes its active target directly from the prepared active chemistry:

```text
Q_active = sum_{i in A} formal_charge_i(prepared molecule)
Q_model = Q_active + sum_F f_F
Q_original = sum_i formal_charge_i(original imported molecule)  # audit provenance only
```

`Q_active` is computed from unselected atoms, without using selected atoms' imported formal charges or
subtracting prescribed charges from the original sum. The prescribed values define fixed sites and their
external field; the modeled molecular total includes them exactly once. `Q_original` remains unchanged
provenance and can differ from `Q_model`. For example, an active neutral protein and a prescribed Mg charge
of +2 have `Q_active = 0` and `Q_model = +2`, whether the selected Mg atom's imported formal charge is 0
or +2. An active prepared protein formal charge of -5 with the same source has `Q_model = -3`.

This is ChargeFW's explicit independent-active-charge model interpretation; the cited embedding equations
support adding a fixed external potential to the active Hamiltonian, but do not validate ChargeFW's active
charge policy or its quantitative accuracy for metal sites. MG and CA CCD records list monatomic formal
charges of +2; a polyatomic CCD component's formal charge is not a partial-charge distribution [16]. The
point-charge embedding literature supports adding the external MM potential to the one-electron/core
Hamiltonian, not inferring the active charge from imported source formal charges [17].

SQE starts from zero and conserves zero charge in each active connected component; SQE+q0 uses prepared
formal-charge seeds, which address SQE's long-zwitterion limitation [15]. SQE+qp normalizes fitted seeds
globally to the active target [1]; subsequent component-local transfers preserve each component's
normalized seed sum, which need not equal its formal-charge total. The peptide paper reports Table 3 BA
typing training-subset RMSDs of 0.0302 for SQE, 0.0188 for SQE+q0, and 0.0133 for SQE+qp. The subset is
part of a 60 di-/tripeptide dataset split 80/20, not all 60 molecules used as training data. Heterogeneous-
set RMSDs favor SQE over q0 (0.0233 vs 0.0290 and 0.0282 vs 0.0350) [1]. These dataset-specific metrics
are not evidence that one variant is generally superior or that any variant is accurate for protein-metal
embedding.
Fixed divalent-ion models can omit induction and charge transfer [18]; this workflow does not add either
response, predict protein protonation, or perform chemical preparation. Input atom formal charges remain
unchanged.

### 4.2 EEM

Current source: [src/methods/builtin/eem.cpp](src/methods/builtin/eem.cpp).

```text
H_ii = B_i
H_ij = kappa / r_ij
chi_i = A_i
phi_i = sum_F kappa * f_F / r_iF
RHS_i = -A_i - phi_i
```

Use the selected parameter set's existing common `kappa`. Do not substitute a universal 14.4 Coulomb
factor. EEM's numeric units/scaling are encoded by that set; for example `EEM_Baek1991` has
`kappa = 0.529177249`. No ion A/B parameter is required for this prescribed source coupling.

The bare point interaction is singular at coincidence. Define validation for nonfinite coordinates,
coincident sources/active atoms, and unphysically close contacts; do not silently clamp distances.

### 4.3 Shared SQE-family solver

Current sources: [sqe.cpp](src/methods/builtin/sqe.cpp),
[sqeq0.cpp](src/methods/builtin/sqeq0.cpp), and [sqeqp.cpp](src/methods/builtin/sqeqp.cpp).
Original SQE and the polypeptide charge-model work provide additional background [14,15]; the linked
implementation determines the exact parameter and sign conventions used here.

Let T be the bond-by-atom incidence matrix, p split charges, s initial charges,
D = diag(H), and K the diagonal bond-hardness matrix. The current implementation is:

```text
(T H T^T + K) p = T [ -chi - H*s + D*s ]
q = s + T^T p
```

With an active-only graph and fixed sources:

```text
(T_A H_AA T_A^T + K) p = T_A [ -chi_A - (H_AA-D_A)*s_A - phi_A ]
q_A = s_A + T_A^T p
```

Preserve the existing `+D*s` seed correction. The code's stored-electronegativity convention gives
`-chi` on the RHS; follow the implementation rather than copying a paper's sign notation blindly.

For an isolated SQE+q0 atom, its incidence column is zero and its seed is its formal charge. The existing
full algebra therefore already keeps it fixed while its off-diagonal seed interactions affect active
transfers. Current parameter classification and pair construction still demand atom parameters. The
new fixed-source path would remove those unnecessary self-parameter requirements and define coupling.

Ordinary SQE still conserves zero charge in each active connected component. Embedding must not silently
make charged active components acceptable. SQE+q0 preserves component formal seed totals; SQE+qp preserves
component sums of its normalized seeds, which need not equal component formal-charge sums.

### 4.4 SQE point-source/active-Gaussian coupling

The current SQE off-diagonal interaction is:

```text
G_ij(r) = erf(r / sqrt(2*w_i^2 + 2*w_j^2)) / r
```

For a point source, set its width to zero:

```text
G_iF(r) = erf(r / (sqrt(2)*abs(w_i))) / r     when w_i != 0
G_iF(r) = 1/r                               when w_i == 0
```

Only the active atom's existing width is needed. A point source does not mean using bare `1/r` against
a Gaussian active atom. A finite-width source would require an explicitly chosen additional width.
Some stored widths are negative; preserve the squared-width behavior.

For nonzero active width, the mathematical zero-distance limit is `sqrt(2/pi)/abs(w_i)`. The current
implementation divides directly by distance and does not implement that limit. Decide explicitly whether
coincident inputs are rejected or an analytic limit is implemented; finite mathematics is not evidence
that overlapping nuclei are scientifically valid.

There is no extra 14.4 factor or EEM common `kappa` in this SQE kernel. SQE bond `kappa` is a different
quantity in the split-charge matrix.

### 4.5 SQE+qp seed normalization

Current normalization reads fitted `q0` for every atom:

```text
s_i = a_i + (Q_total - sum_j a_j) / N
```

An isolated ion therefore stays at its normalized fitted seed, not necessarily its formal charge.
Adding an Mg row with `q0 = 2` would not guarantee Mg2+ in a protein-containing target.

For the proposed model:

```text
s_i = a_i + (Q_A - sum_(j in A) a_j) / N_A     for active atoms
s_F = f_F                                    fixed separately
```

No fixed-atom `q0` is read. Normalize active seeds once for the source active target, including before
fragmentation. Do not switch to component-local normalization: that is a separate scientific change.
Changing a seed budget changes both baseline charges and the `-(H-D)*s` driving term; a post-solve charge
correction cannot repair incorrect normalization.

### 4.6 Method scope and follow-up source inspection

All 22 implemented methods were inspected during the follow-up. The user identified EEM, the SQE family,
and QEq as the main priorities; ABEEM and EQeq are relevant additions, with EQeq+C conditional. The
findings below concern mathematical/implementation feasibility, not tested embedding accuracy. No
embedding solver experiments or code changes were performed.

Freezing a charge and choosing its source kernel are separate decisions. Exact substitution with a
finite screened kernel can still require source hardness. Taking a source-hardness limit is a new,
explicit coupling convention; fixed charge alone does not imply infinite hardness.

| Method | Assessment for this objective |
| --- | --- |
| EEM | Primary; existing `kappa/r` coupling needs no source-specific parameters |
| SQE, SQE+q0, SQE+qp | Primary; explicit point-source/Gaussian coupling and existing seed constraints |
| QEq | Primary; choose an ion-parameter-free source kernel explicitly |
| ABEEM | Strong additional candidate; existing `k/r` acts on atom and bond-charge sites |
| EQeq | Additional candidate; broad elemental coverage already avoids fitted-ion-row gaps |
| EQeq+C | Conditional; settle post-solve correction behavior at the active/fixed boundary |
| Other implemented methods | Outside the requested scope; do not develop new response variants here |

#### ABEEM

In [abeem.cpp](src/methods/builtin/abeem.cpp), nonincident interactions use the selected set's `k/r`.
For isolated fixed ions, subtract `sum_F k*f_F/r_iF` from active atomic RHS entries and
`sum_F k*f_F/r_bF` from active bond-site RHS entries. Bond sites lie at the existing covalent-radius-
weighted centers. No ion atom or bond parameter is required. Driving only atomic sites is incomplete.
The constraint sums atom-site and bond-site charges; output assigns half each bond charge to each end.
Validate source distances to bond centers as well as nuclei. Incident atom-bond coefficients need not
be symmetric, so describe this as restriction of the implemented linear equations, not automatically
as a symmetric energy-Hessian construction. General atomic fragment sources are a chosen representation,
not necessarily exact freezing of an ABEEM fragment's underlying atom and bond sites.

#### QEq

In [qeq.cpp](src/methods/builtin/qeq.cpp), `H_ii = hardness_i` and `RHS_i = -chi_i`. Subtract the
source potential and constrain the active total. All six finite screened kernels require source
hardness, but their formal `hardness_F -> +infinity` limits remove that requirement for positive active
hardness and `r > 0`:

| Overlap option | Limiting active/source kernel |
| --- | --- |
| Nishimoto-Mataga | `14.4/r` |
| Nishimoto-Mataga-Weiss | `17.28/r` |
| Ohno | `14.4/r` |
| Ohno-Klopman | `14.4/sqrt(r*r + (14.4/(2*hardness_i))^2)` |
| DasGupta-Huzinaga (default) | `14.4/r` |
| Louwen-Vogt | `14.4/r` |

These are empirical kernel limits, not proof of a common physical point-charge density interpretation.
In particular, Ohno-Klopman retains active-site softening and Nishimoto-Mataga-Weiss retains its factor
of 1.2. A universal `14.4/r` would not reproduce every limit. The default kernel's limit is a concrete
candidate for parameter-free ion coupling, not an approved final choice. Do not implement the limit by
passing infinite numerical parameters through the existing kernel.

#### EQeq and EQeq+C

[eqeq.cpp](src/methods/builtin/eqeq.cpp) derives electronegativity and hardness from built-in elemental
IP/EA data, with a special hydrogen EA. This provides broad coverage, not guaranteed usable data for
every element or validated coordinated-ion response. Ordinary EQeq still lets the ion charge vary;
elemental coverage does not guarantee Mg stays at +2.

Its constrained solve accepts the same negative potential RHS shift. With
`a = sqrt(hardness_i*hardness_F)/14.4`, the existing pair kernel is:

```text
G_iF(r) = 8.64 * [1/r + exp(-a*a*r*r) * (2*a - a*a*r - 1/r)]
```

Exact finite-screening substitution still needs elemental source hardness. Its infinite-source-hardness
limit at positive distance is `8.64/r`, not `14.4/r`. A different external Coulomb prefactor would define
a different coupling and needs justification.

[eqeqc.cpp](src/methods/builtin/eqeqc.cpp) adds all-pairs geometry corrections after EQeq:

```text
q_i = x_i + sum_(j != i) (Dz_i-Dz_j)*exp[-alpha*(r_ij-R_i-R_j)]
```

Freezing base values `x_F` does not freeze final values `q_F`. Overwriting final source values can break
the total. The simplest candidate is to apply corrections only to active-active pairs after embedded
EQeq; those corrections cancel pairwise in the active total and require no source `Dz`. This deliberately
omits cross-boundary correction terms and must be named as an embedding variant. An exact restriction of
the full correction pipeline instead requires offsetting fixed base values by their corrections and
adjusting the active base budget; then the base field is not generated directly by the prescribed final
charges. No correction-boundary policy was selected.

#### Why the remaining methods are out of scope

SMP/QEq has a suitable iterative equalization structure, but its bundled set lacks ordinary protein
elements. SFKEEM's `2*sqrt(B_i*B_F)/cosh(sigma*r)` requires source `B` and has no useful parameter-free
Coulomb limit. TSEF and DENR admit potential-driven extensions, while PEOE/MPEOE/GDAC and other graph
models would need new coupling assumptions or response interpretations. An unbonded ion has no effect
in the current PEOE-family bond-transfer algorithms. Formal and Dummy have no active response at all.
These observations do not justify expanding this feature into a general method-development project.

## 5. Proposed Application Architecture

### 5.1 Shared request policy, method-specific capability

Prefer a **general calculation/assessment option**, backed by explicit per-method capability and coupling,
over unrelated method options or a fourth execution mode. A name such as `fixed_ion_embedding` is only
illustrative. The internal concept should allow prescribed atomic sources without requiring molecular
fragment presets in the first release.

Embedding is an independent scientific approximation. `full`, `cutoff`, and `cover` remain execution
policies on the active subsystem. A full solve with embedding is not an unmodified full calculation of
the original all-variable-charge system.

Keep toolkit data and preparation policy out of `core::Molecule`. The facade owns the original input,
active input, fixed source data, mappings, and reusable-plan lifetimes. Methods remain stateless and
receive request-specific environment data through their calculation input or a method-scaled potential
derived from it. Raw source data versus precomputed potentials is an implementation choice still to make.

### 5.2 Assessment order

1. Resolve an explicit source selection independently of candidate methods/parameter sets.
2. Validate source identities, charges, coordinates, graph-boundary constraints, and per-target scope.
3. Construct active targets and fixed environments while preserving source indices and conformers.
4. Reject methods that do not support the requested embedding treatment.
5. Evaluate scientific prerequisites on the active target and the embedding-specific requirements.
6. Classify parameters on active atoms and active bonds only, using the unchanged strict matcher.
7. Assess execution availability and select concrete plans normally.
8. Execute the exact assessed partition; reconstruct full source-ordered results once afterward.

An absent Mg parameter should no longer reject an explicitly frozen Mg, but any missing active parameter
must still reject the candidate. Do not automatically freeze atoms merely because a particular candidate
lacks parameters: that would compare different physical systems during automatic selection.

The same policy must reach explicit assessment, direct calculation, automatic selection, reusable plans,
native APIs, Python, and CLI. Bindings/adapters must not duplicate the scientific partitioning policy.

### 5.3 Input and safety contracts to settle

- Default off; no silent changes to existing calculations.
- Initial ion policy: explicit selected atoms, optionally a convenience selection of isolated atoms
  with nonzero supplied formal charges. Isolation is a supplied-graph property, not proof of chemistry.
- Do not infer Mg2+ merely from atomic number or silently remove a selected atom's bonds.
- `core::Atom` stores an integer formal charge defaulting to zero, without a presence flag. Missing and
  explicitly zero charges are not distinguishable there. A neutral fixed atom needs an explicit contract.
- CCD-based ionic charge assignment can be an explicit preparation option using component identity,
  not an unannounced inference. Distinguish FE/FE2 and CU/CU1; validate element/identity correspondence.
- Frozen molecular fragments need a complete finite per-atom charge vector, atom mapping, net-charge
  expectation/tolerance, charge-model provenance, and geometry/conformer policy.
- Reject or explicitly resolve source-total/target-budget inconsistencies. Do not silently renormalize
  source charges or infer a fragment's partial-charge distribution from its net formal charge.
- Keep remaining disconnected active components together in the original target. Splitting them would
  change globally coupled methods and SQE+qp normalization.
- Recheck active-system neutrality requirements. Removing a counterion can expose a charged active target.
- Define all-fixed targets explicitly, including assignment multiplicity and provenance. Do not imply
  that an arbitrary empirical method ran on an empty active graph.
- Scope sources to the intended molecular target; never couple unrelated collection records implicitly.
- Preserve source atom order, molecule/conformer identity, diagnostics, progress mapping, and fixed values.
- Preserve cancellation, thread safety, and lifetime guarantees; avoid global mutable environment state.

### 5.4 Reduced execution

Fixed sources outside an active fragment must still contribute their prescribed field if the embedding
environment is meant to be the same in all execution modes. Precomputing the source potential for each
active atom and projecting it into fragments is one possible design; retain method-specific units/kernels.
For a small number of sources, direct evaluation costs roughly active-atoms times sources.
Resource assessment should account for both the active solve and source-field construction; a large
frozen environment is not computationally free merely because its atoms are no longer solve variables.

Build source reference charges and conserved groups on the active graph. Preserve the original active
SQE+qp normalization through fragments. Apply active conservation corrections before reinserting fixed
charges; sources must never be shifted by those corrections or double-counted as both fragment atoms
and external sources.

Keeping the source field complete does not restore missing active-active interactions or cut bond
transfers. Finite-radius execution remains approximate. Cover chooses source-order pivots, so arbitrary
atom-permutation invariance is not an existing finite-radius cover guarantee.

### 5.5 Provenance

Record the embedding policy independently of execution mode: source atom identities and charge values,
charge provenance, source kernel/width policy, active and original charge budgets, and the frozen-response
assumptions. A structured result field is preferable to a warning alone. An unsupported embedding request
should be a specific rejection, not silently ignored or downgraded to decoupled calculation.

## 6. Code Map for Resumption

| Area | Relevant source |
| --- | --- |
| Request and owned assessment | [assessment.h](include/chargefw/calculation/assessment.h), [assessment.cpp](src/calculation/assessment.cpp) |
| Method input | [calculation_input.h](include/chargefw/methods/calculation_input.h) |
| Method requirements/capabilities | [method_requirements.h](include/chargefw/methods/method_requirements.h) |
| Parameter prerequisites | [parameter_prerequisites.cpp](src/methods/parameter_prerequisites.cpp) |
| Shared execution dispatch | [calculation.cpp](src/calculation/calculation.cpp) |
| Full target budget | [full_execution.cpp](src/calculation/full_execution.cpp) |
| Fragment references/corrections | [reduced_execution.cpp](src/calculation/reduced_execution.cpp) |
| Fragment execution | [cutoff_execution.cpp](src/calculation/cutoff_execution.cpp), [cover_execution.cpp](src/calculation/cover_execution.cpp) |
| Source/target execution mapping | [target_execution.h](src/calculation/target_execution.h) |
| Effective result | [calculation.h](include/chargefw/calculation/calculation.h) |
| Ordered assignments | [charge_collection.h](include/chargefw/charges/charge_collection.h) |
| Gemmi atom import | [structure_import.cpp](src/adapters/gemmi/structure_import.cpp) |
| Gemmi bond policies | [bonds.cpp](src/adapters/gemmi/bonds.cpp) |
| Public classifier semantics | [docs/PARAMETERS.md](docs/PARAMETERS.md) |
| Input preparation boundary | [docs/FORMATS.md](docs/FORMATS.md), [docs/PROJECT.md](docs/PROJECT.md) |
| Existing numerical/fragment tests | [test_sqeqp.cpp](tests/cpp/methods/test_sqeqp.cpp), [test_reduced_execution.cpp](tests/cpp/calculation/test_reduced_execution.cpp) |

Assessment currently has no embedding policy; method input has no external-charge environment. Whole
molecule classification currently happens before calculation. Fixing only the numerical solver will
not make a missing-element target applicable. The reduced executor also contains SQE+qp reference
normalization and must stay consistent with the full path.

## 7. Live PDBe Coverage Research

### 7.1 Definitions and snapshot

Queried live on 2026-10-03; latest indexed release was 2026-09-30. Endpoint:
<https://www.ebi.ac.uk/pdbe/search/pdb/select>.

The index contains **entity documents**, not one document per entry. Aggregated all 609,911
`status:REL` documents, grouped by `pdb_id`, then independently reproduced key counts with server-side
`unique(pdb_id)` facets. Unreleased/withdrawn/obsolete documents must not enter the archive denominator.

"Protein only" below means all indexed polymer entity types are `Protein`, not "contains a protein".
The function filter `number_of_polymers == number_of_protein_chains`, together with `molecule_type:Protein`,
matched that aggregation. `has_modified_residues:N` is a standard-residue metadata proxy, not a complete
coordinate-wide canonical-20 audit.

The water whitelist was explicitly `HOH DOD WAT`. Ordinary water is usually absent from `compound_id`,
but exceptional records include water IDs, especially DOD. These are water-disregarded, not water-free,
cohorts. Production selection should follow its documented water recognition, not assume the search
index and Gemmi selection have identical alias behavior.

### 7.2 Main counts

| Cohort | Unique entries |
| --- | ---: |
| Released PDBe entries | 260,320 |
| Entries containing Protein entities | 255,010 |
| All indexed polymer entities are Protein | 237,677 |
| Above, with no annotated modified residues | 205,882 |
| Baseline: no indexed nonwater compounds | 46,724 |
| At least one NA/K/MG/CA/CL, no other nonwater compounds | 6,380 |
| Baseline plus those five-ion candidates | 53,104 |
| Any of the five ions, other compounds permitted | 51,631 |
| Expanded-list ion(s), no other nonwater compounds | 14,367 |
| Baseline plus expanded-ion candidates | 61,091 |
| Any expanded-list ion, other compounds permitted | 75,321 |

The last seven rows use the 205,882-entry no-modification proxy. The increases over its no-extra baseline
are 13.65% and 30.75%, respectively. These are not measured increases in runnable or accurate calculations.

Expanded list, explicitly selected rather than exhaustive:

```text
NA K MG CA CL ZN MN FE FE2 CU CU1 CO NI CD HG SR CS BA LI RB BR IOD F
```

`IOD`, not `I`, is the CCD iodide ID. All 61,091 combined baseline-plus-expanded entries had
`has_carb_polymer:N` in an independent facet check.

| Ion | Other compounds allowed | Within five-ion-only cohort | Sole nonwater component ID |
| --- | ---: | ---: | ---: |
| CA | 13,450 | 2,398 | 2,060 |
| CL | 16,974 | 1,946 | 1,179 |
| MG | 18,006 | 1,427 | 1,098 |
| NA | 10,370 | 1,312 | 806 |
| K | 3,340 | 343 | 257 |
| Deduplicated union | 51,631 | 6,380 | 5,400 |

Ion rows overlap. Exactly 980 strict entries contain multiple selected ion IDs. Ca-or-Mg occurs in 3,757
strict entries; 68 have both. Multiple copies of a single ion still count as one component ID.

Expanded-list strict per-ion counts, also overlapping:

| ID | Entries | ID | Entries |
| --- | ---: | --- | ---: |
| ZN | 4,837 | IOD | 248 |
| CA | 2,882 | FE2 | 208 |
| CL | 2,576 | HG | 192 |
| MG | 1,846 | CO | 179 |
| NA | 1,555 | BR | 164 |
| MN | 621 | CU1 | 128 |
| FE | 575 | CS | 40 |
| CU | 460 | BA | 32 |
| K | 424 | SR | 15 |
| CD | 341 | RB | 12 |
| NI | 304 | LI | 10 |
| F | 9 | | |

### 7.3 Exclusions and more realistic candidates

Of 51,631 five-ion-containing proxy entries, 45,251 contain additional nonwater compounds. Selected
overlapping occurrences: GOL 7,491; EDO 5,607; SO4 5,533; NAG 3,712; ADP 2,447; ACT 2,397;
PO4 2,308; PEG 2,263; ATP 1,979; GDP 1,684; BMA 1,519; MAN 1,389; HEM 1,261.

These extras are neither automatically removable nor automatically unsupported. Functional ligands may
be essential. Some additives and ligands already match existing parameter sets, as the CCD tests show.

Two later queries deliberately relaxed different restrictions:

- Allow annotated modified component IDs, ignoring the no-modified flag, but no other extras beyond
  the five ions and water: **7,857 candidates**, 1,477 beyond the initial 6,380.
- Retain the no-modified flag but allow `ATP ADP GTP GOL EDO PO4 SO4` as additional components:
  **11,465 candidates**, of which **1,951** contain ATP, ADP, or GTP.

These cohorts overlap and must not be added. Neither proves whole-target parameter coverage. The second
cohort has at least one of the five ions, not necessarily Mg in every nucleotide-containing entry.

For the modified-residue relaxation, group protein-only five-ion documents by entry and evaluate:

```text
compound_ids - {HOH,DOD,WAT} - modified_compound_ids <= {NA,K,MG,CA,CL}
```

The broad input comprised 108,453 entity documents / 60,919 entries. Union component fields across each
entry before applying the test. Among the 1,477 additions, overlapping modified IDs included MSE 653,
ACE 140, SEP 101, NH2 81, CSO 64, CRO 42, TPO 37, LLP 35, PCA 31, PTR 22, CGU 20, CSX 19,
CME 18, KCX 17, MLY 16, NLE 14, DAL 13, DPN 13, and ORN 13. In particular, the substantial MSE
subset remains a selenium-coverage problem for the sets tested.

### 7.4 Examples checked against entry data and coordinates

| Entry | Contents relevant to this question | Strict five-ion candidate |
| --- | --- | --- |
| [1EBH](https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/1ebh) | Enolase, 2 MG, 2 CL, 507 waters | Yes |
| [1EXR](https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/1exr) | Calmodulin, 5 CA, 178 waters | Yes |
| [1CLL](https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/1cll) | Calmodulin, CA, ethanol, water | No: EOH |
| [1ONE](https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/1one) | Enolase, MG, PEP, 2PG, water | No: functional ligands |

1EBH and 1EXR have canonical observed protein residue names and no deposited H atoms. Their ions have
`_atom_site.pdbx_formal_charge = ?`, despite nonzero charges in CCD definitions. Current Gemmi import
copies `atom->charge`; an isolated/nonzero-formal-charge selection alone therefore does not guarantee
that these deposited files will select ions. Explicit chemical preparation/charge assignment is needed.

### 7.5 Reproduce the count queries

This code performs read-only network queries. Responses are live, so later counts may differ.
The complemented regular expression matches component IDs outside the allowed set; its negation excludes
entries with any such component. Explicit `OR` is necessary because the endpoint's default operator is AND.

```python
import json
import requests

URL = "https://www.ebi.ac.uk/pdbe/search/pdb/select"
BASE = "status:REL AND molecule_type:Protein AND has_modified_residues:N"
PROTEIN_ONLY = "{!frange l=0 u=0}sub(number_of_polymers,number_of_protein_chains)"
WATER = {"HOH", "DOD", "WAT"}
FIVE = {"NA", "K", "MG", "CA", "CL"}
EXPANDED = FIVE | set("ZN MN FE FE2 CU CU1 CO NI CD HG SR CS BA LI RB BR IOD F".split())
LIGANDS = set("ATP ADP GTP GOL EDO PO4 SO4".split())

def any_of(ids):
    return " AND compound_id:(" + " OR ".join(sorted(ids)) + ")"

def only(ids):
    return " AND -compound_id:/~(" + "|".join(sorted(ids)) + ")/"

def count(extra=""):
    response = requests.get(URL, params={
        "q": BASE + extra,
        "fq": PROTEIN_ONLY,
        "rows": 0,
        "wt": "json",
        "json.facet": json.dumps({"entries": "unique(pdb_id)"}),
    }, timeout=120)
    response.raise_for_status()
    data = response.json()
    if "error" in data:
        raise RuntimeError(data["error"])
    facets = data["facets"]
    if "entries" in facets:
        return facets["entries"]
    if facets.get("count") == 0:
        return 0
    raise RuntimeError("PDBe response lacks the requested unique-entry count")

baseline = count(only(WATER))
print("proxy", count())                         # 205882
print("baseline", baseline)                    # 46724
for label, ions in [("five", FIVE), ("expanded", EXPANDED)]:
    strict = count(any_of(ions) + only(ions | WATER))
    print(label, strict, baseline + strict, count(any_of(ions)))
    # five: 6380, 53104, 51631; expanded: 14367, 61091, 75321
extra = any_of(FIVE) + only(FIVE | WATER | LIGANDS)
print("with_selected_ligands", count(extra))     # 11465
print("with_nucleotides", count(extra + any_of({"ATP", "ADP", "GTP"})))  # 1951
```

Do not use `response.numFound` as the entry count: the strict five-ion query returned 8,700 documents
but 6,380 entries. For the independent full export, use `q=status:REL`, `rows=50000`,
`sort=entry_entity asc`, and increment `start` by the returned page size. Fields used were:

```text
pdb_id,entry_entity,molecule_type,compound_id,bound_compound_id,modified_compound_id,
has_modified_residues,modified_residue_flag,has_carb_polymer,number_of_polymers,
number_of_protein_chains,release_date
```

`compound_id` includes modified polymer residues and carbohydrate components. `bound_compound_id` can
omit some branched carbohydrate IDs and is not an adequate replacement. The inspected modified-residue
flag `has_modified_residues` was entry-level, whereas `modified_residue_flag` was entity-level. The
stored `molecule_sequence` field did not support an existence query in this endpoint; do not pretend
that a sequence regex certified the canonical-20 condition.

Before explicit water handling, the baseline and strict-five counts were 46,704 and 6,374. DOD-containing
records account for the increases to the final 46,724 and 6,380. Do not reuse the preliminary counts.

### 7.6 Independent RCSB and redundancy checks

RCSB live entry queries returned 255,456 protein-containing entries and 69,271 containing any of the five
ions. `selected_polymer_entity_types = "Protein (only)"` gave 223,469 entries, 56,465 with a selected ion,
and 54,511 with `nonpolymer_entity_count = 0`. These are different filters/index snapshots, not substitutes
for the PDBe no-modification cohorts. Do not combine their denominators.

RCSB sequence-identity grouping supplied a separate broad redundancy check:

| Protein entities from entry cohort | Entity hits | Reported 90% grouping count | Ungrouped entities |
| --- | ---: | ---: | ---: |
| All protein-containing entries | 561,836 | 89,395 | 9,288 |
| Any selected-five-ion entry | 232,989 | 33,062 | 2,654 |
| Protein-only entries | 351,587 | 76,958 | 7,773 |
| Protein-only entries with selected ions | 101,779 | 24,429 | 1,778 |

The search used `return_type: polymer_entity`, a protein-entity predicate, experimental results, and
`group_by: {aggregation_method: sequence_identity, similarity_cutoff: 90}` with representatives and one
returned row. Preserve `total_count`, `group_by_count`, and `ungrouped_count` separately. Entry-level ion
co-presence does not establish that every returned protein entity binds an ion. These groups are not
genes, unique biological proteins, or sequence-deduplicated counts of the strict 6,380 candidates.
The exact RCSB text predicates were `entity_poly.rcsb_entity_polymer_type` equal to `Protein`, optionally
`rcsb_entry_info.selected_polymer_entity_types` equal to `Protein (only)`, and optionally
`rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id` in the five CCD IDs. Use `exact_match`
for the first two and `in` for the third, joined with an `and` group. Request options were
`results_content_type: [experimental]`, `group_by_return_type: representatives`, and
`paginate: {start: 0, rows: 1}` in addition to the grouping specification.
See <https://search.rcsb.org/index.html> and its current request/response schemas before reproducing.

## 8. Biological Significance of the Candidate Pool

The following live PDBe intersections used the strict five-ion query, unique entry counts, and the
specified annotation predicates. Groups overlap; annotations and ion presence alone do not prove an
individual site's physiological role.

| Group | Predicate | Entries | Relevant ion overlap |
| --- | --- | ---: | --- |
| EF-hand domain pair | `interpro_accession:IPR011992` | 500 | 470 CA |
| Cadherin | `pfam_accession:PF00028` | 91 | 90 CA |
| S100 | `pfam_accession:PF01023` | 47 | 43 CA |
| Annexin | `pfam_accession:PF00191` | 27 | 26 CA |
| C-type lectin domain | `pfam_accession:PF00059` | 42 | 32 CA |
| Ion-transport domains | `pfam_accession:(PF00520 OR PF07885)` | 152 | 85 K, 71 CA, 16 NA |
| Enolase activity | `ec_number:4.2.1.11` | 10 | 9 MG |
| PP2C | `interpro_accession:IPR015655` | 20 | 18 MG |
| Inorganic pyrophosphatase | `ec_number:3.6.1.1` | 16 | 9 MG |

A name-defined calmodulin union gave 191 entries, 181 with CA, using exact names `Calmodulin`,
`Calmodulin-1`, `Calmodulin-2`, `Calmodulin-3`, and `Calmodulin-7`. This is not an exhaustive curated family.
GO calcium binding (`GO:0005509`) gave 867 entries / 812 with CA; Mg binding (`GO:0000287`) gave 273 / 196
with MG. There were 2,771 EC-annotated entries. Repetition is substantial: 493 entries had a lysozyme
activity annotation. Missing annotations are not evidence of nonfunctional ions.

Important mechanistic examples:

- **4AQR:** calcium-loaded calmodulin regulates a plasma-membrane calcium pump and calcium homeostasis [10].
- **2WD0:** calcium binding affects cadherin-23 rigidity; a deafness-associated mutation affects affinity [11].
- **5TD9:** human gamma-enolase with MG/CL; reviewed UniProt P09104 identifies Mg's catalytic and
  dimer-stabilizing roles: <https://rest.uniprot.org/uniprotkb/P09104.txt>.
- **2XZV / 2Y09:** PPM/PP2C phosphatase metal-site variants probe the requirement for a third metal [12].
- **2Q6A:** channel calcium binding, permeation, and block involve backbone carbonyl coordination [13].
- **13CM:** a modified-residue candidate contains adenovirus hexon and coagulation factor II with CA/CGU;
  this illustrates a Gla-containing system excluded by a canonical-only filter, not a demonstrated
  ChargeFW calculation: <https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/13cm>.

Use entry `summary`, `molecules`, and `publications` endpoints under
`https://www.ebi.ac.uk/pdbe/api/pdb/entry/` to inspect examples. A counterexample to naive role assignment
is 7SZ8: its search record identifies CA as the ligand of interest and crystallization with 0.1 M MgCl2.
Ion presence is not synonymous with native ion dependence. Calcium-containing enolase 4ROP is another
reason not to infer the native catalytic cofactor from family and deposited ion alone.

## 9. Local CCD and Native Parameter-Coverage Experiments

### 9.1 Input and method

Local input: repository-root `components.gz`, supplied by the user and untracked during research.
Do not add this large file to a commit merely because this document references it.

SHA-256:

```text
05becffae2d3de5b9a2840832f8415eb9f9a672cb8a5170c52141c63d409c174
```

Streamed 51,346 blocks through Gemmi 0.7.4, one block at a time. Do not load the entire decompressed
dictionary unnecessarily. A replacement CCD release can change results; the hash identifies this input.

Used existing Python bindings from `build/gcc-debug/python`. System Python could import Gemmi, but
needed that `PYTHONPATH` to import ChargeFW. Tests called native assessment with explicit method/set,
`parameter_matching="strict"`, and `execution="full"`. Distinct **nonphysical placeholder coordinates**
only satisfied the coordinate prerequisite; no solver or charge-accuracy measurement was performed.

| Label | Exact parameter-set ID |
| --- | --- |
| PUB | `SQEqp_Schindler2021_PUB_pept` |
| CCD | `SQEqp_Schindler2021_CCD_gen` |
| DTP | `SQEqp_Schindler2021_DTP_small` |
| NEEMP | `EEM_NEEMP_ccd2016_npa` |
| Cheminf | `EEM_Cheminf_b3lyp_npa` |
| Baek | `EEM_Baek1991` |

PUB uses neighbor-element `bonded` atom types and typed bonds. CCD/DTP use highest-bond-order `hbo`
atom types plus typed bonds. Residue names are not classification keys. Formal charge is not directly
part of these atom keys, but changing protonation changes the graph and can change the types.
Permissive matching lowers selected bond-order classifications; it does not invent elements or neighbors.

### 9.2 Linked-residue results

Constructed `Gly-X-Gly` with NH3+ / COO- outer termini. For X, checked named
`N-CA-C(=O)-OXT` connectivity, replaced one N-H with an upstream peptide bond, and replaced C-OXT(H)
with a downstream peptide bond. Other CCD H atoms, bond orders, and formal charges were preserved unless
an explicit variant is identified. ALA/GLY controls passed all sets.

| Central X / state | PUB | CCD | DTP | NEEMP | Cheminf | Baek |
| --- | --- | --- | --- | --- | --- | --- |
| NLE, DAL, DPN | Pass | Pass | Pass | Pass | Pass | Pass |
| ORN, neutral side-chain NH2 | Fail | Pass | Pass | Pass | Pass | Pass |
| ORN, explicitly protonated side-chain NH3+ | Pass | Pass | Pass | Pass | Pass | Pass |
| HYP, CGU, KCX, MLY, unchanged CCD side chains | Fail | Pass | Pass | Pass | Pass | Pass |
| CGU, both side-chain carboxylates deprotonated | Fail | Pass | Pass | Pass | Pass | Pass |
| CME, CSO | Fail | Pass | Pass | Pass | Pass | Fail |
| SEP, TPO, PTR, LLP | Fail | Pass | Fail | Pass | Pass | Pass |
| CSX | Fail | Pass | Fail | Pass | Pass | Fail |
| MSE | Fail | Fail | Fail | Fail | Fail | Fail |

These are complete parameter-assessment successes/failures for the constructed graphs, not residue-name
guarantees in arbitrary proteins. Not all phosphate-containing residue protonation states were tested.

Selected PUB failures, with `C[CCCH]` meaning an atom with those immediate neighbor elements:

- Neutral ORN: single `C[CHHN]-N[CHH]`; adding a side-chain N-H and assigning N +1 rescues coverage.
- HYP: single `C[CCHH]-C[CCHO]` and `C[CCHO]-C[CHHN]`; atom types exist but bond pairs do not.
- CGU: two single `C[CCCH]-C[COO]` bonds remain missing even after deprotonating both carboxyls.
- CME: single `C[CHHS]-C[CHHO]`.
- CSO: `S[CO]`, `O[HS]`, and associated bonds.
- CSX: `S[CHO]`, `O[S]`, and associated bonds.
- KCX: `C[NOO]` and associated bonds.
- MLY: methyl `C[HHHN]` and associated bonds.
- MSE: Se and Se-adjacent environments. It cannot be treated as an isolated fixed ion inside the residue.
- SEP/TPO/PTR/LLP: phosphorus environments, among other PUB gaps.

17/20 canonical central residues passed PUB with unchanged CCD side-chain forms. ASP/GLU passed after
explicit carboxyl deprotonation. HIS retained the missing double bond `N[CCH]-C[HNN]` for its tested
CCD tautomer/Kekule representation. Do not turn that result into a blanket assertion that all histidine
representations fail. Disulfide types/bonds are present in PUB; do not exclude disulfides categorically.

PCA passed a hypothetical internal N-acylated `Gly-PCA-Gly` graph. Ordinary biological pyroglutamate is
usually N-terminal; that result does not certify terminal PCA coverage and is not a headline success.

### 9.3 Standalone and cap traps

PUB failed every tested standalone selected CCD component, including canonical amino acids. The isolated
neutral terminal groups differ from internal peptide environments. This does not mean linked residues fail.
CCD passed the selected standalone components except MSE. DTP additionally failed SEP/TPO/PTR/LLP/CSX;
the broader EEM sets passed the selected components except MSE. Baek also failed selected sulfur cases.

An explicit `Ac-X-NH2` construction failed PUB even for ALA because the cap introduces the missing single
bond `C[CNO]-C[CHHH]`. A capping artifact must not be reported as an intrinsic residue failure. Matched
canonical controls are essential before interpreting a modified-residue test.

### 9.4 Ligands and explicit nucleotide anions

| Standalone graph | PUB | CCD | DTP | NEEMP | Cheminf | Baek |
| --- | --- | --- | --- | --- | --- | --- |
| Raw CCD ATP, ADP, GTP | Fail | Pass | Fail | Pass | Pass | Pass |
| Explicit ATP4-, ADP3-, GTP4- | Fail | Pass | Fail | Pass | Pass | Pass |
| GOL, EDO | Fail | Pass | Pass | Pass | Pass | Pass |
| PO4(-3) | Fail | Pass | Fail | Pass | Pass | Pass |
| SO4(-2) | Fail | Pass | Fail | Pass | Pass | Fail |

Raw local CCD nucleotide records have net formal charge zero and protonated phosphate acids. For the
anion tests, removed only H atoms attached to O singly bonded to phosphate P, set each resulting O to -1,
and preserved all other retained atoms and bond orders. Removed four H for ATP/GTP (`HOG2`, `HOG3`,
`HOB2`, `HOA2`) and three for ADP (`HOB2`, `HOB3`, `HOA2`); verified totals -4, -4, and -3.

These are explicitly chosen protonation states, not automatic pH preparation or a claim that every
bound nucleotide adopts that state. They establish that existing broader sets can cover important
ligand chemistry; they do not certify a full metal-protein-ligand assembly or response accuracy.
Do not mix independently fitted PUB/CCD parameter rows to create an unvalidated hybrid set.

### 9.5 Ion inventory

Found 81 charged monatomic CCD component IDs, corresponding to 80 distinct stored element/charge
combinations. ZN/ZN2 duplicate Zn2+. The inventory includes rare/radioactive elements and deuteron D8U;
it is not an inventory of 81 scientifically validated embedding applications.

| Charge | Common CCD IDs |
| --- | --- |
| +1 | NA, K, LI, RB, CS, CU1 |
| +2 | MG, CA, ZN, MN, FE2, CU, CO, NI, CD, HG, SR, BA |
| +3 | FE |
| -1 | CL, BR, IOD, F |

None of these 81 isolated-ion records passed the three tested SQEqp sets or the NEEMP/Cheminf sets.
Baek matched only AL. Organic element coverage does not guarantee isolated-ion type coverage: `hbo`
also changes when an atom has no bonds.

### 9.6 Durable reproduction recipe

The essential experiment can be recreated without the temporary scripts:

1. Stream gzip text, splitting only at CIF `data_` block boundaries; parse each block with Gemmi.
2. Read `_chem_comp_atom.atom_id`, `type_symbol`, `charge`, and `pdbx_leaving_atom_flag`.
3. Read `_chem_comp_bond.atom_id_1`, `atom_id_2`, and `value_order`. The selected tests used explicit
   SING/DOUB/TRIP mapped to 1/2/3; reject unfamiliar orders rather than guessing.
4. Preserve all CCD H atoms, charges, and bond orders for standalone tests.
5. Construct linked graphs explicitly as described above, with charged terminal Gly controls. Do not
   assume that dropping every atom flagged as leaving is sufficient to build a polymer.
6. For side-chain variants, remove the relevant O-H and set O -1 for carboxylates; for ORN add one NE-H
   and set NE +1. For nucleotide variants, modify phosphate O-H only.
7. Use `chargefw.Molecule`, then native strict/full assessment with an explicit method and parameter set.
8. Record `bool(assessment.plans)` and all rejection messages, including source atom/bond indices.
9. Never run charge solvers with the nonphysical assessment-only placeholder coordinates.

Minimal standalone/inventory reproduction core (run with the existing built/installed bindings; see
[DEVELOPMENT.md](DEVELOPMENT.md) for supported build setup):

```python
import gzip
import gemmi
import chargefw

def ccd_blocks(path):
    with gzip.open(path, "rt") as stream:
        lines = []
        for line in stream:
            if line.startswith("data_") and lines:
                yield gemmi.cif.read_string("".join(lines)).sole_block()
                lines = []
            lines.append(line)
        if lines:
            yield gemmi.cif.read_string("".join(lines)).sole_block()

def molecule_from_block(block):
    atoms = block.get_mmcif_category("_chem_comp_atom.")
    bonds = block.get_mmcif_category("_chem_comp_bond.")
    names = atoms["atom_id"]
    index = {name: i for i, name in enumerate(names)}
    orders = {"SING": 1, "DOUB": 2, "TRIP": 3}
    return chargefw.Molecule(
        atomic_numbers=[gemmi.Element(x).atomic_number for x in atoms["type_symbol"]],
        formal_charges=[int(x) for x in atoms["charge"]],
        atom_names=names,
        bonds=[(index[a], index[b], orders[o]) for a, b, o in zip(
            bonds.get("atom_id_1", []), bonds.get("atom_id_2", []),
            bonds.get("value_order", []))],
        # Assessment only: these are not physical molecular geometries.
        coordinates=[[float(i), float(i % 3), float(i % 7)] for i in range(len(names))],
        name=block.name,
    )

selected = {"NLE", "HYP", "MSE", "CGU", "ATP", "ADP", "GTP", "GOL", "EDO", "PO4", "SO4"}
ion_ids = []
block_count = 0
for block in ccd_blocks("components.gz"):
    block_count += 1
    atoms = block.get_mmcif_category("_chem_comp_atom.")
    if len(atoms.get("atom_id", [])) == 1 and int(atoms["charge"][0]) != 0:
        ion_ids.append((block.name, atoms["type_symbol"][0], int(atoms["charge"][0])))
    if block.name not in selected:
        continue
    molecule = molecule_from_block(block)
    for suffix in ("PUB_pept", "CCD_gen", "DTP_small"):
        pid = "SQEqp_Schindler2021_" + suffix
        assessment = chargefw.assess(
            molecule, method="sqeqp", parameter_set=pid,
            parameter_matching="strict", execution="full",
        )
        print(block.name, pid, bool(assessment.plans))
        for rejection in assessment.rejections:
            for issue in rejection.issues:
                print(" ", issue.message)
print("CCD blocks", block_count, "charged monatomic components", len(ion_ids))
```

This minimal snippet intentionally does not claim that standalone results reproduce linked-residue
results. The exact graph construction and protonation choices above are needed for that experiment.
The durable recipe reproduces chemical graphs up to atom ordering and their coverage outcomes, not the
original numeric rejection indices. Preserve explicit constructed graphs when index-exact comparisons
are required; equivalent choices of deleted N-H names or appended-atom order can change those indices.

Temporary research artifacts were left in `/tmp/kilo/`:

```text
ccd_coverage.py                 ccd_selected.json
ccd_coverage_results.json       ccd_variants.py
ccd_variants_results.json       ccd_ligands.py
ccd_ligand_results.json         ccd_ions.json
ccd_ion_coverage.py             ccd_ion_coverage_results.json
ccd_nucleotide_anions.py        ccd_nucleotide_anion_results.json
```

The JSON records include explicit constructed graphs and native rejection messages. These files are
ephemeral and not part of the repository; the findings and reconstruction rules in this document must
remain sufficient if they disappear. The main script used six sets, including the three EEM sets above.

## 10. Validation and Implementation Sequence

### 10.1 Focused initial scope

1. Approve a concrete public request/capability contract and explicit source-charge policy.
2. Prototype full execution for EEM and the shared SQE family; include QEq after approving its
   ion-parameter-free source kernel. ABEEM and EQeq are additional candidates, not prerequisites.
3. Add owned active/frozen preparation, active-only classification, diagnostics, and source reassembly.
4. Validate mathematical invariants and a small scientific benchmark before broad claims.
5. Add cutoff/cover propagation and validate separately. Until supported, reject such requested combinations.
6. Expose identical semantics through native, Python, and CLI, with structured provenance.
7. Treat molecular frozen sources and curated presets as optional follow-up work, not prerequisites.

Do not make PEOE support, all metal models, covalent fragment cuts, automatic protonation, or new fitted
metal parameters prerequisites for the initial feature. Do not silently alter published parameter sets.

### 10.2 Numerical and behavioral tests

- Fixed values remain exact; active and combined totals follow the requested budget.
- Missing fixed-ion parameters cease to matter, while missing active parameters still fail.
- With embedding disabled, calculations reproduce existing behavior; unsupported embedding methods/policies
  reject clearly.
- Compare a small explicitly constrained quadratic system with the active-plus-potential formulation.
- Vary unused source self parameters in a reference construction: active results must not depend on them.
- A stable symmetric two-site active fixture responds with more negative induced charge on the side
  nearer a positive source. Do not demand this monotonicity for every atom in arbitrary molecules.
- Source reversal and superposition follow the linear model at fixed active geometry and budget.
- A uniform potential does not redistribute charge within a constrained group; far-field differential
  response vanishes appropriately.
- Verify method-specific scaling, signed widths, zero widths, coincidence policy, and invalid coordinates.
- Check SQE+qp active-only seed normalization and preserve `-H*s + D*s` exactly.
- Test isolated/multiple/negative ions, no active atoms, disconnected active components, and source-total
  conflicts. Test complete molecular-source vectors and boundaries if that interface is included.
- Verify atom permutations and bond endpoint reversals for full solves and all mapping paths.
- Verify multiple molecules/conformers, source order, source lifetime, plan reuse, cancellation, and threads.
- In reduced execution, include distant sources outside fragment radii, avoid double counting, keep
  sources out of conservation corrections, and recover full results at full-system radius.
- For finite-radius cover, test mapping and convergence without inventing a permutation-invariance promise.

### 10.3 Scientific benchmark

Use representative small complexes/site models rather than starting with whole-protein QM:

- Hydrated Mg2+/Ca2+ and carboxylate/amide coordination, with deliberate water/protonation treatment.
- Calcium-binding EF-hand examples and a magnesium phosphatase-site model.
- Mg-nucleotide models using an explicit ligand protonation state and a parameter set that actually matches it.
- If molecular frozen sources are included, sulfate/phosphate tests with justified atomic source charges.
- Monovalent-ion and distant-ion controls.

Separate two validation questions: (a) does EEM/SQE reproduce the response to a prescribed point-charge
field; (b) how well does that prescribed-field approximation represent the real ion-containing system?
Compare embedded and ion-omitted calculations on the same prepared active geometry. Assess metal-site
and first-shell atom errors separately from whole-protein averages. Include distance scans, conditioning,
overpolarization, sensitivity to source-charge model, electrostatic potentials/dipoles when relevant,
and held-out coordination environments. Conserved charge alone is not accuracy evidence.

Use a reference charge scheme consistent with the intended observable. NPA agreement and electrostatic
potential reproduction are different objectives. For reaction chemistry, strong covalency, redox, or
proton transfer, the frozen model may be outside its useful scope even if classification passes.

To measure actual software coverage gain, prepare a representative PDB corpus with one documented
import/protonation/water policy and assess it before/after embedding. Separate failures due to ions,
active parameter types, missing atoms/H, geometry, resource policy, and numerical execution. Record both
entry counts and a meaningful nonredundant grouping; metadata composition counts are only a screening step.

Follow the repository validation cadence in [DEVELOPMENT.md](DEVELOPMENT.md) and [AGENTS.md](AGENTS.md):
focused gcc-debug first; coherent full debug suite; clang cross-compiler checks; release checks for
numerics; sequential conservative ASan/UBSan for mappings/views/indexing; clang-tidy for implementation.
No builds or charge solver tests were run as part of this preliminary research; native assessment used
an existing build. Rebuild and establish a clean baseline before implementing.

## 11. Open Decisions and Recommendation

The follow-up confirmed the objective and priorities: prescribe ion charges without fitted ion response
parameters, retain the protein's response, prioritize EEM/SQE-family/QEq, and avoid overengineering.
ABEEM and EQeq remain relevant; EQeq+C is conditional. The user explicitly accepts an ion-only feature
if molecular fragments add substantial complexity. No commitment to a fragment library or automatic
partial-charge assignment was made.

The recommended minimal design uses mapped atomic sources and the original target's conformer positions,
with no active/fixed graph bonds. This keeps a later explicit molecular charge vector inexpensive without
requiring its public interface now. Component identification, source-charge assignment, and field
coupling are separate responsibilities; do not build a generic source-provider framework solely for
possible future presets. Parameter-covered small molecules can remain active.

Unresolved decisions include the public policy name/shape, source selection syntax, explicit CCD ionic
charge assignment, QEq source-kernel choice, whether explicit molecular vectors are exposed initially,
coincidence handling, all-fixed output semantics, potential ownership, and result provenance. If EQeq+C
is included, its correction boundary also needs approval. Prefer in-target sources initially; independent
external charge clouds and fragment-charge generation are deferred recommendations, not required scope.
Public API and automatic chemistry choices still need approval before implementation.

The recommendation is a focused, explicit embedding capability, not new fitted metal parameters or
external-field variants of every implemented method.
The scientific motivation is not merely "avoid an error for Mg": it is retaining the electrostatic
environment of functionally important ions while reusing already covered protein, modified-residue,
and ligand chemistry. For calcium regulation and Mg-nucleotide workflows this could be an enabling
capability. Its accuracy must be established for selected use cases, not inferred from either PDB
counts or the fact that related approximations are standard elsewhere.

## References

1. Schindler et al. (2021), *Optimized SQE atomic charges for peptides accessible via a web application*.
   DOI: <https://doi.org/10.1186/s13321-021-00528-w>. Full text:
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC8243439/>. Source parameterization and seed normalization.
2. Senn and Thiel (2009), *QM/MM methods for biomolecular systems*.
   DOI: <https://doi.org/10.1002/anie.200802019>. Established active/environment modeling framework.
3. Zinovjev (2023), *Electrostatic Embedding of Machine Learning Potentials*.
   DOI: <https://doi.org/10.1021/acs.jctc.2c00914>. Full text:
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC10061678/>. Explicit embedding discussion; its response model
   is not a validation of bare EEM/SQEqp.
4. Cisneros, Piquemal, and Darden (2006), *Quantum Mechanics/Molecular Mechanics Electrostatic Embedding
   with Continuous and Discrete Functions*. DOI: <https://doi.org/10.1021/jp062768x>. Full text:
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC2656107/>. Short-range representation/polarization effects,
   including Mg2+-water examples.
5. Tzeliou, Mermigki, and Tzeli (2022), *Review on the QM/MM Methodologies and Their Application to
   Metalloproteins*. DOI: <https://doi.org/10.3390/molecules27092660>. Full text:
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC9105939/>. Electrostatic versus polarizable embedding.
6. Rick, Stuart, and Berne (1994), *Dynamical fluctuating charge force fields: Application to liquid water*.
   DOI: <https://doi.org/10.1063/1.468398>. Environment-responsive, electronegativity-based charges in MD.
7. Joung and Cheatham (2008), *Determination of alkali and halide monovalent ion parameters for use in
   explicitly solvated biomolecular simulations*. DOI: <https://doi.org/10.1021/jp8001614>. Full text:
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC2652252/>. Coulomb/Lennard-Jones ion models and water-specific
   calibration, not a parameter-free universal ion model.
8. Li and Merz (2017), *Metal Ion Modeling Using Classical Mechanics*.
   DOI: <https://doi.org/10.1021/acs.chemrev.6b00440>. Full text:
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC5312828/>. Model choices, limitations, polarization, and CT.
9. Kashefolgheta and Vila Verde (2017), *Developing force fields when experimental data is sparse:
   AMBER/GAFF-compatible parameters for inorganic and alkyl oxoanions*.
   DOI: <https://doi.org/10.1039/c7cp02557b>. Sulfate/phosphate modeling precedent; abstract verified,
   publisher full text was inaccessible during this research. No detailed charge-generation prescription
   was inferred from the abstract.
10. *A bimodular mechanism of calcium control in eukaryotes* (2012).
    DOI: <https://doi.org/10.1038/nature11539>. Calmodulin/calcium-pump regulation; PDB 4AQR.
11. Sotomayor et al. (2010), *Structural determinants of cadherin-23 function in hearing and deafness*.
    DOI: <https://doi.org/10.1016/j.neuron.2010.03.028>. PDB 2WD0 and calcium-dependent mechanics.
12. *A third metal is required for catalytic activity of the signal-transducing protein phosphatase M
    tPphA* (2011). DOI: <https://doi.org/10.1074/jbc.M109.036467>. PDB 2XZV/2Y09.
13. *Structural insight into Ca2+ specificity in tetrameric cation channels* (2007).
    DOI: <https://doi.org/10.1073/pnas.0707324104>. PDB 2Q6A.
14. Nistor et al. (2006), *A generalization of the charge equilibration method for nonmetallic materials*.
    DOI: <https://doi.org/10.1063/1.2346671>. Original SQE formulation.
15. Verstraelen et al. (2012), *Assessment of atomic charge models for gas-phase computations on
    polypeptides*. DOI: <https://doi.org/10.1021/ct200512e>. SQE-family polypeptide/initial-charge context.
16. RCSB Protein Data Bank Chemical Component Dictionary records: magnesium `MG` formal charge +2,
    <https://files.rcsb.org/ligands/view/MG.cif>; calcium `CA` formal charge +2,
    <https://files.rcsb.org/ligands/view/CA.cif>. Monatomic component records, not a polyatomic
    partial-charge assignment.
17. External molecular-mechanics point-charge potential in the one-electron/core Hamiltonian.
    DOI: <https://doi.org/10.1002/jcc.27532>. This supports external-field coupling, not ChargeFW's
    active-charge budget or metal-site accuracy.
18. Fixed divalent-ion model limitations regarding induction and charge transfer.
    DOI: <https://doi.org/10.1021/ct400751u>;
    <https://pmc.ncbi.nlm.nih.gov/articles/PMC3960013/>.

Primary archive sources: <https://www.ebi.ac.uk/pdbe/search/pdb/select>,
<https://www.ebi.ac.uk/pdbe/api/>, and <https://search.rcsb.org/>.
Counts, local CCD coverage, native-assessment results, and proposed scientific validity are separate
levels of evidence throughout this document.
