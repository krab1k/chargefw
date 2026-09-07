# Python recipes

These executable examples demonstrate complete Python workflows. The [Python package guide](../PYTHON.md)
contains shorter API snippets, while the [format reference](../FORMATS.md) defines shared molecular input
and charge output behavior.

## Inspect imported molecules

[`inspect_molecules.py`](inspect_molecules.py) reads any supported molecular format and reports the source
identity, dimensions, conformers, and total formal charge of every imported record. Structural selection,
bond, and conformer policies are explicit options.

```bash
python docs/recipes/inspect_molecules.py input.sdf --format sdf
```

## Calculate an SDF collection

[`calculate_sdf_collection.py`](calculate_sdf_collection.py) reads a multimolecule SDF, assesses one
explicit QEq policy for the complete collection, calculates source-mapped charges, and writes result JSON
with requested and effective provenance. One inapplicable record makes the policy inapplicable to the
collection; the recipe does not silently choose a separate model per record.

```bash
python docs/recipes/calculate_sdf_collection.py input.sdf result.json
```

Python molecular output is generated from normalized molecules and does not preserve SDF data fields.
Result JSON is used here because it retains the complete calculation record.

## Charge a Gemmi document

[`charge_gemmi_document.py`](charge_gemmi_document.py) attaches calculated charges to a caller-owned Gemmi
mmCIF document. Original mmCIF input is enriched without reconstructing the document; PDB input is
necessarily converted to a new mmCIF document. The input must already have appropriate hydrogens, formal
charges, coordinates, and explicit or template-supported topology.

```bash
pip install "chargefw[gemmi]"
python docs/recipes/charge_gemmi_document.py input.cif charged.cif --format mmcif
```

The recipe uses full QEq with `QEq_original` by default. Full execution can be expensive for large
structures, and the recipe reports rather than overrides resource warnings.

## Analyze conformer-dependent charges

[`analyze_rdkit_conformer_charges.py`](analyze_rdkit_conformer_charges.py) generates an unoptimized,
deterministic RDKit ETKDG ensemble, calculates QEq charges, and ranks atoms by charge standard deviation
over the equally weighted generated conformers.

```bash
pip install "chargefw[rdkit]"
python docs/recipes/analyze_rdkit_conformer_charges.py "CCO" --conformers 20
```

This is a sensitivity analysis over one generated coordinate set, not a thermodynamic ensemble or a
claim that the ranked atoms are chemically significant hotspots. The useful analysis function also
accepts an externally prepared RDKit conformer ensemble.

## Compare bundled parameter sets

[`compare_parameter_sets.py`](compare_parameter_sets.py) assesses every bundled parameter set for one
method and one molecular geometry, reports inapplicable sets, and compares executable source-aligned
charge vectors against one reference.

```bash
python docs/recipes/compare_parameter_sets.py molecule.sdf --format sdf --method eem
```

Differences demonstrate parameter sensitivity only. Parameter sets may target different molecular
domains or reference charge schemes; the recipe does not rank their scientific accuracy.
