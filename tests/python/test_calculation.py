"""Focused calculation facade and result-model checks."""

import numpy as np

import chargefw


def water(conformers=1):
    coordinates = np.array(
        [
            [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
            [[0.0, 0.0, 0.1], [0.96, 0.0, 0.1], [-0.24, 0.93, 0.1]],
        ][:conformers]
    )
    return chargefw.Molecule(
        [8, 1, 1],
        bonds=np.array([[0, 1, 1], [0, 2, 1]], dtype=np.int32),
        coordinates=coordinates[0] if conformers == 1 else coordinates,
        source_name="water.sdf",
        record_index=2,
        record_id="water-record",
        atom_ids=["O", "H1", "H2"],
        conformer_ids=["model-a"] if conformers == 1 else ["model-a", "model-b"],
    )


calculator = chargefw.Calculator()
options = chargefw.CalculationOptions(method="eem", execution="full")
assessment = calculator.assess(water(), options)
assert assessment.executable
assert assessment.execution_policy == {
    "mode": "full",
    "radius": None,
    "charge_correction": "none",
}
assert assessment.report["applicability"] == assessment.applicability
selected_index = assessment.applicability["selected_candidate_index"]
selected_candidate = assessment.applicability["applicable"][selected_index]
assert selected_candidate["method_id"] == "eem"

result = assessment.calculate()
assert result.status == "success"
assert len(result.assignments) == 1
assignment = result.assignments[0]
assert assignment.molecule_index == 0
assert assignment.conformer_index == 0
assert assignment.source.record_id == "water-record"
assert assignment.source_atom_ids == ("O", "H1", "H2")
assert assignment.source_conformer_id == "model-a"
assert assignment.values.dtype == np.float64
assert assignment.values.flags.c_contiguous
assert not assignment.values.flags.writeable
assert np.isclose(assignment.values.sum(), 0.0)
assert result.requested["method_id"] == "eem"
assert result.effective["method_id"] == "eem"
assert result.effective["parameter_set_id"] == selected_candidate["parameter_set_id"]
assert result.effective["execution_policy"]["mode"] == "full"
assert result.timings["applicability_seconds"] >= 0.0
assert result.timings["computation_seconds"] >= 0.0
try:
    options.method_options["eem"] = {"unexpected": True}
except TypeError:
    pass
else:
    raise AssertionError("calculation options must be immutable")

geometry_independent = calculator.calculate(
    chargefw.Molecule([8, 1, 1]),
    chargefw.CalculationOptions(method="formal", execution="full"),
)
assert geometry_independent.status == "success"
assert len(geometry_independent.assignments) == 1
assert geometry_independent.assignments[0].conformer_index is None

try:
    assessment.calculate()
except RuntimeError:
    pass
else:
    raise AssertionError("assessment must be one-shot")

multi_result = calculator.calculate(
    chargefw.MoleculeCollection([water(2)]),
    chargefw.CalculationOptions(method="qeq", execution="full"),
)
assert multi_result.status == "success"
assert [item.conformer_index for item in multi_result.assignments] == [0, 1]
assert [item.source_conformer_id for item in multi_result.assignments] == ["model-a", "model-b"]

reduced = calculator.calculate(
    water(),
    chargefw.CalculationOptions(
        method="eem", execution="cutoff", radius=8.0, charge_correction="none"
    ),
)
assert reduced.status == "success"
assert reduced.effective["execution_policy"] == {
    "mode": "cutoff",
    "radius": 8.0,
    "charge_correction": "none",
}

covered = calculator.calculate(
    water(), chargefw.CalculationOptions(method="eem", execution="cover", radius=8.0)
)
assert covered.status == "success"
assert covered.effective["execution_policy"]["mode"] == "cover"

no_plan = calculator.assess(
    chargefw.Molecule([8, 1, 1], bonds=[[0, 1, 1], [0, 2, 1]]),
    chargefw.CalculationOptions(method="qeq", execution="full"),
).calculate()
assert no_plan.status == "no_executable_plan"
assert no_plan.assignments == ()
try:
    no_plan.raise_for_status()
except chargefw.NoExecutablePlanError as error:
    assert error.result is no_plan
else:
    raise AssertionError("no-plan result must raise its typed exception")

try:
    chargefw.CalculationOptions(execution="full", charge_correction="uniform")
except ValueError:
    pass
else:
    raise AssertionError("full execution must reject charge correction")
