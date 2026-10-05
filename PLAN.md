# Fixed Ions Plan

## Objective

Provide opt-in fixed ions for EEM, SQE, SQE+q0, and SQE+qp. Explicitly selected,
fixed-ion components create a field to which active atoms respond; their charges are not appended
to an unperturbed calculation. The existing solver, cutoff, and cover machinery is reused.

The priority is SQE, SQE+q0, and SQE+qp. EEM remains supported as a secondary method, not as a reason
to prioritize broader legacy parameter coverage. Fixed ions do not require fitted parameters for
the selected fixed components; ordinary strict parameter matching remains in force for the active
graph. Molecular components such as sulfate remain active.

Fixed ions default off. Python `assess()` and `calculate()` accept named monatomic ions independently of
candidate method or parameter coverage:

```python
calculate(molecule, fixed_ions=["MG", "CA"])
```

The Python request snapshots names into an immutable tuple. `None` and an empty sequence disable
resolution. `COMMON_IONS` and `ALL_IONS` are immutable tuples selecting the common subset and complete
bundled monatomic-ion catalog. Template charges are authoritative; no imported-first charge model is used.

## Decisions

- Resolve named component instances to original atom indices in the native preparation/facade layer.
  Bindings and applications pass intent; they do not implement scientific policy.
- Do not expose `MonatomicComponentCharge` or require users to construct records containing element,
  charge, and template vectors. The approved Python API uses named ions; internal indexed source values
  remain numerical plumbing. Effective result provenance retains its existing fixed-charge representation.
- Keep input atom formal charges unchanged. The active target charge is the formal-charge sum of
  unselected atoms in the prepared active molecule; fixed charges contribute exactly once to the
  modeled total. Keep the active target internal; results retain only indexed fixed values. See [METALS.md](METALS.md)
  and the owning [project contract](docs/PROJECT.md) for rationale and implementation boundaries.
- Fixed component templates are bundled, immutable, versioned data. A named ion template assigns its
  formal charge, validates component identity, element, and complete atom coverage, and applies only
  to a matching instance in the same target. PDB/mmCIF use must not require network access or depend
  on a component-supplied charge table being present. Different identities such as FE and FE2 are distinct.
- Do not infer oxidation state from element alone. Do not silently select components because a
  candidate lacks their parameters. Selection is caller intent and is independent of candidate
  coverage.
- Recognize selected monatomic components by catalog ID, exactly one atom of the expected element,
  and absence of incident graph bonds. Atom-name matching is unnecessary. Do not infer ions from
  graph singletons or guess charge from element alone.
- Keep molecular components active. Fixed-ion provenance identifies one indexed atom per instance;
  generic imported component membership continues to represent complete molecular components.
- Expose the same immutable native resolution through native assessment, Python plans, and CLI requests.
  Keep resolved source indices and values explicit in results.
- Keep execution policy distinct from the fixed-charge approximation. Existing full, cutoff, and cover
  machinery supplies fixed sources to every relevant solve, including sources beyond a fragment radius.

## Named-Ion Resolution

Structural imports own explicit selected-residue component membership in import metadata. The private
named-ion resolver consumes that partition rather than rebuilding instances from hierarchy labels; original
labels remain audit metadata. The partition and complete import metadata remain owned with native records
and Python source mappings. Its private catalog and common-ion subset are defined in the owning
[project contract](docs/PROJECT.md#assessment-and-execution). Coordination connectivity remains limited
to bonds represented in the imported core graph.
Template and fixed-ion matching share label-first canonical names so alternate author names do not change
component identity, while source mappings continue to preserve both namespaces.

The user-facing request selects named ions. Native indexed charges and method-level positional point
sources remain general numerical inputs, independent of the adapter catalog.

## Workflow

Astra owns design and review; Luna owns substantial implementation, with Astra handling small local edits.
Follow the design-review and qualification workflow in [AGENTS.md](AGENTS.md); agent review is not user
approval. Commands belong in [DEVELOPMENT.md](DEVELOPMENT.md), implementation history in git, and only
current decisions and unfinished work here.
