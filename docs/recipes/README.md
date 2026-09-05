# Python recipes

These executable examples are both user documentation and API design checks. Keep them direct, concise,
and runnable; awkward recipe code is a reason to reconsider the public API rather than hide it behind
recipe-specific helpers.

## Conformational charge hotspots

[`conformational_charge_hotspots.py`](conformational_charge_hotspots.py) generates an RDKit molecule from
SMILES, adds explicit hydrogens, embeds conformers, calculates QEq charges, and reports the atoms with the
lowest and highest charge variance across conformers.

```bash
pip install "chargefw[rdkit]"
python docs/recipes/conformational_charge_hotspots.py "CCO" --conformers 20
```

The recipe intentionally produces no serialized molecule or charge output.

## Charge a Gemmi structure

[`charge_structure_with_gemmi.py`](charge_structure_with_gemmi.py) loads a PDB or mmCIF structure with
Gemmi, converts it to ChargeFW molecules, calculates charges, attaches them to a caller-owned Gemmi mmCIF
document, and writes the document.

```bash
pip install "chargefw[gemmi]"
python docs/recipes/charge_structure_with_gemmi.py input.cif charged.cif
```
