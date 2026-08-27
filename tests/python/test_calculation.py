"""Focused calculation facade and result-model checks."""

import numpy as np
from pathlib import Path

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
options = chargefw.CalculationOptions(
    method="eem", execution=chargefw.ExecutionSelectionKind.FULL
)
assessment = calculator.assess(water(), options)
assert assessment.executable
assert assessment.execution_policy == chargefw.ExecutionPolicy(
    mode=chargefw.ExecutionMode.FULL,
    radius=None,
    charge_correction=chargefw.ChargeCorrectionPolicy.NONE,
)
assert assessment.report.applicability == assessment.applicability
selected_index = assessment.applicability.selected_candidate_index
selected_candidate = assessment.applicability.applicable[selected_index]
assert selected_candidate.method_id == "eem"
assert selected_candidate.execution_assessments[0].mode is chargefw.ExecutionMode.FULL

result = assessment.calculate()
assert result.status is chargefw.ExecutionStatus.SUCCESS
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
assert assignment.values.base is not None
assert np.isclose(assignment.values.sum(), 0.0)
assert result.requested is options
assert result.requested.method == "eem"
assert result.effective.method_id == "eem"
assert result.effective.parameter_set_id == selected_candidate.parameter_set_id
assert result.effective.execution_policy.mode is chargefw.ExecutionMode.FULL
assert result.timings.applicability_seconds >= 0.0
assert result.timings.computation_seconds >= 0.0
try:
    options.method_options["eem"] = {"unexpected": True}
except TypeError:
    pass
else:
    raise AssertionError("calculation options must be immutable")

geometry_independent = calculator.calculate(
    chargefw.Molecule([8, 1, 1]),
    chargefw.CalculationOptions(
        method="formal", execution=chargefw.ExecutionSelectionKind.FULL
    ),
)
assert geometry_independent.status is chargefw.ExecutionStatus.SUCCESS
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
    chargefw.CalculationOptions(method="qeq", execution=chargefw.ExecutionSelectionKind.FULL),
)
assert multi_result.status is chargefw.ExecutionStatus.SUCCESS
assert [item.conformer_index for item in multi_result.assignments] == [0, 1]
assert [item.source_conformer_id for item in multi_result.assignments] == ["model-a", "model-b"]

reduced = calculator.calculate(
    water(),
    chargefw.CalculationOptions(
        method="eem",
        execution=chargefw.ExecutionSelectionKind.CUTOFF,
        radius=8.0,
        charge_correction=chargefw.ChargeCorrectionPolicy.NONE,
    ),
)
assert reduced.status is chargefw.ExecutionStatus.SUCCESS
assert reduced.effective.execution_policy == chargefw.ExecutionPolicy(
    mode=chargefw.ExecutionMode.CUTOFF,
    radius=8.0,
    charge_correction=chargefw.ChargeCorrectionPolicy.NONE,
)

covered = calculator.calculate(
    water(),
    chargefw.CalculationOptions(
        method="eem", execution=chargefw.ExecutionSelectionKind.COVER, radius=8.0
    ),
)
assert covered.status is chargefw.ExecutionStatus.SUCCESS
assert covered.effective.execution_policy.mode is chargefw.ExecutionMode.COVER

no_plan = calculator.assess(
    chargefw.Molecule([8, 1, 1], bonds=[[0, 1, 1], [0, 2, 1]]),
    chargefw.CalculationOptions(method="qeq", execution=chargefw.ExecutionSelectionKind.FULL),
).calculate()
assert no_plan.status is chargefw.ExecutionStatus.NO_EXECUTABLE_PLAN
assert no_plan.assignments == ()
assert no_plan.applicability.rejected
assert isinstance(
    no_plan.applicability.rejected[0].issues[0].kind,
    chargefw.PrerequisiteIssueKind,
)
try:
    no_plan.raise_for_status()
except chargefw.NoExecutablePlanError as error:
    assert error.result is no_plan
else:
    raise AssertionError("no-plan result must raise its typed exception")

try:
    chargefw.CalculationOptions(
        execution=chargefw.ExecutionSelectionKind.FULL,
        charge_correction=chargefw.ChargeCorrectionPolicy.UNIFORM,
    )
except ValueError:
    pass
else:
    raise AssertionError("full execution must reject charge correction")

try:
    chargefw.CalculationOptions(execution="full")
except TypeError:
    pass
else:
    raise AssertionError("execution strings must not be accepted")

try:
    chargefw.CalculationOptions(charge_correction="none")
except TypeError:
    pass
else:
    raise AssertionError("charge correction strings must not be accepted")

try:
    chargefw.CalculationOptions(max_threads=True)
except TypeError:
    pass
else:
    raise AssertionError("boolean thread limits must not be accepted")

try:
    chargefw.CalculationOptions(cutoff_atom_threshold=False)
except TypeError:
    pass
else:
    raise AssertionError("boolean atom thresholds must not be accepted")

try:
    chargefw.CalculationOptions(cover_atom_threshold=np.iinfo(np.uintp).max + 1)
except ValueError:
    pass
else:
    raise AssertionError("atom thresholds must fit the native size type")

methods = calculator.methods
eem = next(method for method in methods if method.id == "eem")
assert eem.requires_coordinates
assert eem.supports_cutoff
assert eem.supports_cover
assert isinstance(eem.options, tuple)
assert all(option.type is chargefw.MethodOptionType.INTEGER for option in eem.options)
assert chargefw.method_descriptors() == methods

parameter_directory = Path(chargefw.__file__).parent / "_data" / "parameters"
eem_path = next(parameter_directory.glob("EEM_Ouy2009.json"))
parameter_set = chargefw.load_parameter_set(eem_path)
assert parameter_set.descriptor.method_id == "eem"
assert parameter_set.id == parameter_set.descriptor.id
loaded_sets = chargefw.load_parameter_sets(parameter_directory)
assert parameter_set.id in {value.id for value in loaded_sets}
explicit_calculator = chargefw.Calculator([parameter_set])
assert explicit_calculator.parameter_sets == (parameter_set.descriptor,)
explicit_result = explicit_calculator.calculate(water(), options)
assert explicit_result.status is chargefw.ExecutionStatus.SUCCESS
assert explicit_result.effective.parameter_set_id == parameter_set.id
try:
    chargefw.Calculator([])
except ValueError:
    pass
else:
    raise AssertionError("an explicit parameter catalog must not be empty")
