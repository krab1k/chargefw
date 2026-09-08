# Python recipes

These executable examples demonstrate complete Python workflows. Ordinary calculation recipes call
`chargefw.calculate()` directly. The parameter-set comparison recipe uses explicit assessment as an
advanced workflow for inspecting and executing several prepared plans. The
[Python package guide](../PYTHON.md) contains shorter API snippets, while the
[format reference](../FORMATS.md) defines shared molecular input and charge output behavior.

## Calculate a molecular file

[`calculate_file.py`](calculate_file.py) is the shortest complete file-to-charges workflow. It reads any
supported molecular format, calculates one explicitly selected method and parameter set in full, and
prints each source-ordered charge array.

```bash
python docs/recipes/calculate_file.py ethanol.sdf \
    --format sdf \
    --method qeq \
    --parameter-set QEq_original
```

## Inspect imported molecules

[`inspect_molecules.py`](inspect_molecules.py) reads any supported molecular format and reports the source
identity, dimensions, conformers, and total formal charge of every imported record. Structural selection,
bond, and conformer policies are explicit options.

```bash
python docs/recipes/inspect_molecules.py input.sdf --format sdf
```

## Calculate an SDF collection

[`calculate_sdf_collection.py`](calculate_sdf_collection.py) reads a multimolecule SDF, directly
calculates one explicit QEq policy for the complete collection, and writes result JSON with requested and
effective provenance. One inapplicable record makes the policy inapplicable to the collection; the recipe
does not silently choose a separate model per record.

```bash
python docs/recipes/calculate_sdf_collection.py input.sdf result.json
```

Python molecular output is generated from normalized molecules and does not preserve SDF data fields.
Result JSON is used here because it retains the complete calculation record.

For independent automatic policies, use an ordinary per-record loop instead. Each call performs a separate
assessment, so the selected method, parameter set, or execution mode can differ between records:

```python
molecules = chargefw.io.read("input.sdf", format="sdf")

for molecule in molecules:
    result = chargefw.calculate(molecule)
    for assignment in result.assignments:
        print(molecule.name, result.plan.method.id, assignment.values)
```

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

[`compare_parameter_sets.py`](compare_parameter_sets.py) performs one explicit assessment for a method
and molecular geometry, reports rejected parameter sets, and executes each prepared plan to compare its
source-aligned charge vector against one reference. This advanced workflow avoids repeating molecule
preparation and parameter classification for every executable parameter set.

```bash
python docs/recipes/compare_parameter_sets.py molecule.sdf --format sdf --method eem
```

Differences demonstrate parameter sensitivity only. Parameter sets may target different molecular
domains or reference charge schemes; the recipe does not rank their scientific accuracy.
