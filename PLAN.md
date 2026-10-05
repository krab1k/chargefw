# Fixed-Charge Groups Plan

## Objective

Provide opt-in fixed-charge groups for EEM, SQE, SQE+q0, and SQE+qp. Explicitly selected,
fixed-charge components create a field to which active atoms respond; their charges are not appended
to an unperturbed calculation. The existing solver, cutoff, and cover machinery is reused.

The priority is SQE, SQE+q0, and SQE+qp. EEM remains supported as a secondary method, not as a reason
to prioritize broader legacy parameter coverage. Fixed-charge groups do not require fitted parameters for
the selected fixed components; ordinary strict parameter matching remains in force for the active
graph. A parameter-covered component such as sulfate remains active unless the caller selects it.

Fixed-charge groups default off. Python `assess()` and `calculate()` accept named component groups independently of
candidate method or parameter coverage:

```python
calculate(molecule, fixed_charge_groups=["MG", "CA"])
```

The Python request snapshots names into an immutable tuple. `None` and an empty sequence disable
resolution. `COMMON_IONS` and a charge-model selector are not part of the API.

## Decisions

- Resolve named component instances to original atom indices in the native preparation/facade layer.
  Bindings and applications pass intent; they do not implement scientific policy.
- Do not expose `MonatomicComponentCharge` or require users to construct records containing element,
  charge, and template vectors. The approved Python API uses named groups; internal indexed source values
  remain numerical plumbing. Effective result provenance retains its existing fixed-charge representation.
- Keep input atom formal charges unchanged. The active target charge is the formal-charge sum of
  unselected atoms in the prepared active molecule; fixed charges contribute exactly once to the
  modeled total. Preserve the original formal-charge sum as audit provenance. See [METALS.md](METALS.md)
  and the owning [project contract](docs/PROJECT.md) for rationale and implementation boundaries.
- Fixed component templates are bundled, immutable, versioned data. A named ion template assigns its
  formal charge, validates component identity, element, and complete atom coverage, and applies only
  to a matching instance in the same target. PDB/mmCIF use must not require network access or depend
  on a component-supplied charge table being present. Different identities such as FE and FE2 are distinct.
- Do not infer oxidation state from element alone. Do not silently select components because a
  candidate lacks their parameters. Selection is caller intent and is independent of candidate
  coverage.
- First scope: opt-in whole graph-isolated components. Internal bonds are allowed when a molecular
  template is supported; no crossing bond deletion, capping, or active/fixed bond cutting. Isolation
  is only as reliable as the imported graph. Coordination-bond import semantics must be explicit.
- Input formal charges are not a molecular partial-charge distribution. Ions are required now;
  small molecular presets can follow when their charge distributions are justified. HEM is desirable
  but not a first-release requirement. No HEM or other molecular charge distribution is adopted.
- A later explicitly named formal-charge model or scientifically sourced partial-charge model is
  possible, but the exact model and option are unsettled. Do not offer a vague "more chemical" mode,
  arbitrary unknown-charge generator, or preset for every ligand.
- Expose the same immutable native resolution through assess/calculate and reusable Python plans. The CLI
  integration remains future work. Keep provenance explicit: preset ID/version, resolved values, and any
  approximation warning.
- Keep execution policy distinct from the fixed-charge approximation. Existing full, cutoff, and cover
  machinery supplies fixed sources to every relevant solve, including sources beyond a fragment radius.

## Named-Ion Resolution

Structural imports own explicit selected-residue component membership in import metadata. The private
named-ion resolver consumes that partition rather than rebuilding instances from hierarchy labels; original
labels remain audit metadata. The partition and complete import metadata remain owned with native records
and Python source mappings. Its initial private catalog is `NA`, `K`, `MG`, `CA`, `CL`, `ZN`, `FE`, and
`FE2`. Coordination connectivity remains limited to bonds represented in the imported core graph.
Template and fixed-ion matching share label-first canonical names so alternate author names do not change
component identity, while source mappings continue to preserve both namespaces.

The desired user-facing request stays simple: named components plus any explicitly necessary charge
model choice. Public element/charge/template-vector records are rejected. Internal per-atom templates
and indexed values may still feed numerical source coupling.

## Remaining Work

- [ ] Integrate named fixed-charge groups into the CLI without duplicating native resolution policy.
- [ ] Qualify SQE-family priority against independent indexed/manual references. Cover source identity,
  same-target mapping, conformers, original atom order, invalid/incomplete mappings, fixed values,
  active charge totals, and EEM behavior as secondary coverage.
- [ ] Scope multiatom templates separately. Verify complete atom identity/element matching, versioned
  provenance, internal versus crossing bonds, and justified charge distributions before adding presets.
- [ ] Review scientific limits and publish implemented behavior in the owning docs only after the
  feature and its tests are complete.

## Workflow

Astra owns design and review; Luna owns substantial implementation, with Astra handling small local edits.
Follow the design-review and qualification workflow in [AGENTS.md](AGENTS.md); agent review is not user
approval. Commands belong in [DEVELOPMENT.md](DEVELOPMENT.md), implementation history in git, and only
current decisions and unfinished work here.
