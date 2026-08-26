# ChargeFW Delivery Roadmap

This file contains unfinished work only. Implemented state belongs in [PROJECT.md](PROJECT.md), usage
in [README.md](README.md), and implementation rules in [AGENTS.md](AGENTS.md). Remove an item when its
full acceptance criteria are met; do not retain checked history here.

Work is ordered by current product dependency: establish scientific reference evidence, make the native
package consumable and continuously checked, then add Python distribution. ACC III replacement has
separate adoption gates. Optional integrations and distribution formats should be added only for a
concrete supported workflow.

## 1. Scientific validation and ChargeFW2 comparison

The implementation and execution contracts have focused regression coverage. The remaining work is
scientific evidence and release readiness. ChargeFW2 is a previous implementation and useful comparison
material, not an oracle: a difference is a finding to investigate, not evidence by itself that ChargeFW
is wrong. Cutoff and cover remain explicit approximations until their method-specific validation is
complete.

- [ ] Build a reusable ChargeFW2 comparison harness that runs both implementations with identical
  topology, formal charges, conformers, method options, and parameter data. Keep small comparison
  inputs, source-order mappings, observed outputs, and method-specific tolerances in version control.
- [ ] Use the harness to examine every built-in method and every bundled parameter set with a ChargeFW2
  counterpart on at least one compatible molecule. Add neutral, ionic, disconnected, and
  multi-conformer cases where they exercise relevant method behavior. For each difference, investigate
  the equations, parameters, input interpretation, numerical solver, and legacy implementation before
  classifying it as a ChargeFW defect, a ChargeFW2 defect, an intentional change, or unresolved.
- [ ] Investigate ChargeFW2 cutoff behavior for the overlapping EEM/QEq-like methods (`eem`, `eqeq`,
  `eqeqc`, and `qeq`) at representative radii. Compare fragment selection, target charge, correction,
  and failure behavior, but do not treat reproducing the old result as the acceptance criterion. Record
  the scientific rationale for the selected current behavior and any unresolved differences.
- [ ] Establish a representative full-versus-cutoff/cover corpus for the eight reduced-capable methods
  across multiple radii, charge states, disconnected systems, conformers, and at least one practically
  large structure. Report mean/max charge error, charge conservation, runtime, and peak memory, then
  state method-specific accuracy guidance without implying exact equivalence.
- [ ] Review each built-in method's equations, default options, prerequisites, diagnostics, and citation
  against its publication and the archived implementation while adding the numerical comparisons. Use
  the publication and current scientific reasoning to resolve conflicts; the archived implementation is
  evidence, not the final authority.
- [ ] Add targeted robustness cases for the actual dense-solver failure modes that remain untested,
  especially singular or nearly singular systems. Require either finite charge-conserving output or an
  actionable failure with molecule/conformer context; do not build a combinatorial edge-case matrix.

## 2. Native installation and automation

- [ ] Add a lean CI baseline for the supported development platform: formatting, GCC debug/release
  tests, Clang ASan/UBSan, and the native install/consumer smoke test. Add further compiler, platform,
  or package jobs when they protect a supported distribution rather than maintaining an unused matrix.
- [ ] Give bundled parameter data a machine-readable release identifier and retain it in calculation
  provenance. Before 1.0, document how software, result-schema, and parameter-data versions affect
  compatibility and reproducibility.

## 3. Python and toolkit integration

- [ ] Define a synchronous toolkit-neutral Python API accepting atomic numbers, formal charges,
  indexed bonds, source names/identities, and zero or more coordinate arrays; return NumPy charge
  arrays plus mappings, provenance, and diagnostics.
- [ ] Add nanobind bindings over the owned facade and test ownership, exception translation, array
  validation, collection ordering, source mappings, and all-conformer results.
- [ ] Build a base wheel with bundled parameter data, package-resource discovery that needs no
  environment variable, and clean-install tests for a declared initial CPython/platform set. Expand
  the matrix only to platforms the project intends to support.
- [ ] Add a pure-Python `rdkit.Chem.Mol` converter and an explicit charge-attachment helper behind a lazy
  optional dependency. Preserve atom indices, formal charges, supported bonds, and selected conformers;
  perform no implicit sanitization, hydrogen changes, protonation, embedding, or optimization, and do
  not overwrite existing properties unless requested.

## 4. ACC III adoption gates

- [ ] Establish a reproducible ACC III comparison corpus with licensed/reference inputs, parameter
  versions, options, mappings, and method-specific tolerances. Treat licensed ACC III/ChargeFW2 outputs
  as comparison observations, not automatically correct answers.
- [ ] Investigate each material ChargeFW2/ACC III divergence and classify it as a ChargeFW defect, a
  legacy/reference defect, an intentional change, or unresolved. Retain the rationale and only promote
  behavior to a regression requirement when the scientific or product reason is understood.
- [ ] Before replacing the ACC III backend, define acceptance criteria for its intended workflows and
  obtain independent scientific review of the method comparisons and any cutoff/cover accuracy claims.
