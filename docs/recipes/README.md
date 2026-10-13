# Python recipes

These executable examples demonstrate complete Python workflows. Ordinary calculation recipes call
`chargefw.calculate()` directly. The parameter-set comparison recipe uses explicit assessment as an
advanced workflow for inspecting and executing several prepared plans. The
[Python package guide](../PYTHON.md) contains shorter API snippets, while the
[format reference](../FORMATS.md) defines shared molecular input and charge output behavior.

## Calculate a molecular file

[`calculate_file.py`](calculate_file.py) is the shortest complete file-to-charges workflow. It reads any
supported molecular format, calculates one explicitly selected method, and prints one row per atom in
source order with its molecule and conformer indices, atom ID, atom name, element, and charge. A `-`
conformer identifies charges from a geometry-independent method. `--result-json` also writes complete
result JSON with requested and effective provenance.

```bash
python docs/recipes/calculate_file.py molecules.sdf --format sdf --method eem --result-json result.json
```

One method, parameter set, and execution mode apply to the whole file, so one inapplicable record makes
the request inapplicable; the recipe does not silently choose a separate model per record. Its EEM
invocation lets ChargeFW select a bundled parameter set applicable to every record. The
[Python collection guide](../PYTHON.md#reading-molecular-files) explains when to calculate a collection or
its records independently, and the [format reference](../FORMATS.md) defines import and result JSON
behavior.

## Inspect imported molecules

[`inspect_molecules.py`](inspect_molecules.py) reads any supported molecular format and reports the source
identity, dimensions, conformers, and total formal charge of every imported record. Structural selection,
bond, and conformer policies are explicit options. With `--classifier plain`, `hbo`, or `bonded`, it also
counts each record's atoms by element and [parameter classifier](../PARAMETERS.md#classifiers) type,
showing how finely that classifier separates the atoms, for example before grouping calculated charges.

```bash
python docs/recipes/inspect_molecules.py input.sdf --format sdf --classifier bonded
```

## Charge a Gemmi document

[`charge_gemmi_document.py`](charge_gemmi_document.py) attaches charges to the unchanged live Gemmi mmCIF
document used for import. PDB input is first represented as a caller-owned Gemmi mmCIF document. See the
[Gemmi integration contract](../PYTHON.md#toolkit-integrations) for validation and mutation behavior.

```bash
pip install "chargefw[gemmi]"
python docs/recipes/charge_gemmi_document.py input.cif charged.cif --format mmcif
```

The recipe uses EEM by default and lets ChargeFW select an applicable bundled parameter set.

## Calculate with a fixed ion

[`calculate_with_fixed_ions.py`](calculate_with_fixed_ions.py) reads a PDB or mmCIF structure, selects one
or more monatomic ion components as fixed sources, calculates the remaining active atoms, and prints the
resolved ion indices and charges. You can also write the complete result JSON.

```bash
python docs/recipes/calculate_with_fixed_ions.py structure.cif \
    --format mmcif --ion MG --ion CA --result-json charges.json
```

The selected components should be present in the input. The recipe uses `sqeqp` with
`SQEqp_Schindler2021_CCD_gen` by default; override the method or parameter set when the active structure
requires another model.

## Analyze conformer-dependent charges

[`analyze_rdkit_conformer_charges.py`](analyze_rdkit_conformer_charges.py) generates an unoptimized,
deterministic RDKit ETKDG ensemble, calculates EEM charges, and ranks atoms by charge standard deviation
over the equally weighted generated conformers.

```bash
pip install "chargefw[rdkit]"
python docs/recipes/analyze_rdkit_conformer_charges.py "CCO" --conformers 20
```

This is a sensitivity analysis over one generated coordinate set, not a thermodynamic ensemble or a
claim that the ranked atoms are chemically significant hotspots. The useful analysis function also
accepts an externally prepared RDKit conformer ensemble.

## Write charges to RDKit molecules

[`charge_rdkit_sdf.py`](charge_rdkit_sdf.py) reads an SDF with RDKit, calculates all records with one
method, parameter set, and execution mode, attaches the charges as the `ChargeFWPartialCharge` atom
property, and records the model as molecule properties. The output SDF stores the charges as an
`atom.dprop.ChargeFWPartialCharge` property list, which RDKit restores as atom properties when it reads
the file. Records that RDKit cannot read are reported and skipped. Every remaining record must convert,
so records need explicit hydrogens.

```bash
pip install "chargefw[rdkit]"
python docs/recipes/charge_rdkit_sdf.py input.sdf charged.sdf --method eem
```

Because the converted records form one collection, ChargeFW selects a parameter set applicable to all of
them, and a record the selected model cannot handle makes the calculation fail rather than switching
models. Pass `--parameter-set` to choose the model explicitly. The recipe kekulizes aromatic rings so
ChargeFW receives explicit single and double bonds; see the
[RDKit integration contract](../PYTHON.md#toolkit-integrations).

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
