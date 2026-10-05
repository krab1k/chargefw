# Fixed Ions: Rationale and Evidence

This note preserves compact evidence for [PLAN.md](PLAN.md), not an implemented-behavior contract.
Follow [AGENTS.md](AGENTS.md); build procedures belong in [DEVELOPMENT.md](DEVELOPMENT.md). User-facing
architecture and charge accounting contracts belong in [docs/PROJECT.md](docs/PROJECT.md), serialization
in [docs/FORMATS.md](docs/FORMATS.md), and interface behavior in the native, Python, or CLI document.

## Decision Context

Explicitly selected fixed ions let active atoms respond to a prescribed electrostatic field
without requiring fitted metal parameters. They do not automatically freeze every component lacking
candidate parameters: sulfate and other molecular components remain active. Priority is
SQE/SQE+q0/SQE+qp; EEM remains secondary, not a justification for prioritizing old parameter coverage.
Existing full, cutoff, and cover machinery supplies the fixed sources to relevant solves.

This is a one-way fixed ion treatment: active charges respond; source polarization and
active/fixed charge transfer are absent. QM/MM point-charge methods provide an analogy, not validation
of ChargeFW response. Organic response parameters are not certified for metal fields; quantitative
metal-site accuracy is unvalidated.

Charge accounting uses the independent active-charge model:

```text
Q_active = sum(formal charges of unselected atoms in prepared active molecule)
Q_model = Q_active + sum(prescribed fixed charges)
```

Input formal charges remain unchanged. Preserve mappings and conformers and insert fixed values exactly
in original atom order. An input formal-charge sum differing from the modeled total is not a
reason to renormalize active atoms or sources. Only the internal active target is needed for calculation.

Bundled immutable named-component templates provide deterministic charges when structure input omits
them. Validate component identity, element identity, and complete atom mapping; do not infer oxidation
state from element alone (FE and FE2 are distinct). PDB/mmCIF use must not require a network lookup or
component-supplied input charge table. A multiatom component's formal total does not determine its
partial-charge distribution.

## Structure-Index Evidence

### PDBe Snapshot (2026-10-05)

PDBe `/search/pdb/select`, query `status:REL`, reports 260,320 released entries. The `compound_id`
facet counts below are unique entries containing each component ID; counts overlap and do not count
copies or imply runnable calculations. There is no protein or canonical-residue filter.

| Component ID | Entries | Component ID | Entries |
| --- | ---: | --- | ---: |
| MG | 28,940 | ZN | 27,985 |
| CA | 16,623 | SO4 | 28,307 |
| NAG | 16,163 | MSE | 10,295 |
| HEM | 6,473 | ADP | 5,018 |
| ATP | 3,802 | FAD | 3,343 |
| FES | 1,549 | HEC | 1,321 |

The `compound_id` index is incomplete: it misses SF4 even for released entry 10EG, while RCSB reports
SF4. Index absence or a zero facet count is not evidence of true absence. Separate RCSB experimental
entry counts are SF4 2,576, F3S 374, and B12 178. These use a different service/query and have no
claimed common denominator or released-entry filter; do not combine with PDBe.

### Reproduction Shapes

PDBe component distribution (all released entries, no protein filter):

```text
GET https://www.ebi.ac.uk/pdbe/search/pdb/select
  ?q=status:REL
  &rows=0
  &json.facet={"entries":"unique(pdb_id)","components":{"type":"terms",
    "field":"compound_id","limit":150,"sort":"entries desc",
    "facet":{"entries":"unique(pdb_id)"}}}
  &wt=json
```

RCSB count query uses `POST https://search.rcsb.org/rcsbsearch/v2/query`, attribute
`rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id`, operator `in`,
`return_type: entry`, `results_content_type: [experimental]`; read `total_count`:

```json
{
  "query": {
    "type": "terminal",
    "service": "text",
    "parameters": {
      "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
      "operator": "in",
      "value": ["SF4"]
    }
  },
  "return_type": "entry",
  "request_options": {
    "results_content_type": ["experimental"],
    "paginate": {"start": 0, "rows": 5}
  }
}
```

Counts are the 2026-10-05 snapshot; live indices can change.

## CCD And Native Assessment

The user-supplied, untracked `components.gz` had SHA-256
`05becffae2d3de5b9a2840832f8415eb9f9a672cb8a5170c52141c63d409c174`. Gemmi 0.7.4 streamed 51,346
blocks. Existing `build/gcc-debug` Python bindings ran strict, full-execution native assessment. Seven
CCD graphs were checked: HEM, HEC, SF4, FES, F3S, B12, and CNC. Placeholder coordinates served only
assessment; no solver or accuracy measurement was run. CCD records include explicit hydrogens, formal
charges, and bonds; unknown charges defaulted to zero in B12/CNC assessments, limiting their value.

| CCD graph assessment | EEM sets accepting | SQE-family sets accepting |
| --- | ---: | ---: |
| HEM, HEC, SF4, FES, F3S (each) | 2/28 | 0/9 |
| B12, CNC (each) | 0/28 | 0/9 |

The EEM acceptors for the first five are `EEM_Svob2007_cmet2` and `EEM_Svob2007_hm2`. Ordinary SQE
may also fail its neutral-component prerequisite, so its failures are not solely matcher failures.
These are component-graph assessments, not full protein/ligand assembly coverage.

Separately, the tested parameter inventory had potential element coverage for P in 11/28 EEM sets and
3/9 SQE-family sets; Fe in 2/28 EEM and 0/9 SQE sets; Mg in 0/28 and 0/9; and Se or Co in 0/28 and
0/9. These are element-presence counts, not graph matching results. Even ordinary Fe acceptance by
some EEM sets does not establish coverage of the protein's active graph, its bonds, or the complete
parameter set needed by a whole calculation, and it says nothing about charge accuracy.

Broader checked sets accepted standalone sulfate, phosphate, nucleotides, sugars, and additives, but
this was not an all-sets survey or quantified whole-assembly coverage. MSE is commonly covalently linked
in proteins, not an isolated source. The evidence favors removing fixed metal parameter requirements
for the SQE priority while retaining strict matching for active chemistry.

## Charge Records And Boundaries

CCD `MG` has `_chem_comp_atom.charge = +2`. CCD `SO4` assigns S 0, two O atoms 0, and two O atoms -1.
This localized formal-charge pattern is not a justified physical partial-charge distribution; a
symmetric S 0 / four O -0.5 heuristic changes the local field too and is not adopted.

For 1EBH, deposited `_chem_comp_atom` has no charge column and `_atom_site` lacks ionic formal charge.
PDB has no chemical-component dictionary. Bundled templates therefore provide a reliable offline
identity/charge source when inputs omit this information.

Heme instances may have covalent attachments, especially heme c, or imported metal-coordination links;
boundary eligibility depends on the actual graph and import policy. HEC is generally covalently
attached. SF4/FES/F3S cluster identity, ligands, and redox state matter. Graph isolation is only as
reliable as imported connectivity. Ion templates are the initial need; small molecular presets require
justified charge distributions. HEM is desirable but not a first-release requirement, and no heme
charge distribution is adopted. A formal-charge template and a named partial-charge model are distinct
possibilities; their identity and public option remain unsettled.

## References

- Schindler et al. (2021), SQE peptide parameterization and scope. DOI:
  <https://doi.org/10.1186/s13321-021-00528-w>.
- Senn and Thiel (2009), QM/MM methods for biomolecular systems. DOI:
  <https://doi.org/10.1002/anie.200802019>.
- Cisneros, Piquemal, and Darden (2006), external electrostatic potentials and short-range limitations. DOI:
  <https://doi.org/10.1021/jp062768x>.
- Li and Merz (2017), metal-ion modeling, polarization, and charge-transfer limitations. DOI:
  <https://doi.org/10.1021/acs.chemrev.6b00440>.
- RCSB CCD: [MG](https://files.rcsb.org/ligands/view/MG.cif),
  [SO4](https://files.rcsb.org/ligands/view/SO4.cif), and [1EBH](https://files.rcsb.org/view/1EBH.cif).
- PDBe search: <https://www.ebi.ac.uk/pdbe/search/pdb/select>; RCSB search: <https://search.rcsb.org/>.

These sources establish an external-field precedent, not mutual polarization, charge transfer, or
quantitative ChargeFW accuracy for metal-site response.
