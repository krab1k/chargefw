"""Focused calculation facade and result-model checks."""

from pathlib import Path
import unittest
from typing import Any, cast

import numpy as np

import chargefw


def water(conformers: int = 1) -> chargefw.Molecule:
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


class CalculationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.calculator = chargefw.Calculator()
        self.full_eem = chargefw.CalculationOptions(
            method="eem", execution=chargefw.ExecutionSelectionKind.FULL
        )

    def test_assessment_and_full_calculation(self) -> None:
        assessment = self.calculator.assess(water(), self.full_eem)
        self.assertTrue(assessment.executable)
        self.assertEqual(
            assessment.execution_policy,
            chargefw.ExecutionPolicy(
                mode=chargefw.ExecutionMode.FULL,
                radius=None,
                charge_correction=chargefw.ChargeCorrectionPolicy.NONE,
            ),
        )
        self.assertEqual(assessment.report.applicability, assessment.applicability)
        selected_index = assessment.applicability.selected_candidate_index
        if selected_index is None:
            self.fail("executable assessment must select a candidate")
        selected_candidate = assessment.applicability.applicable[selected_index]
        self.assertEqual(selected_candidate.method_id, "eem")
        self.assertIs(
            selected_candidate.execution_assessments[0].mode,
            chargefw.ExecutionMode.FULL,
        )

        result = assessment.calculate()
        self.assertIs(result.status, chargefw.ExecutionStatus.SUCCESS)
        self.assertEqual(len(result.assignments), 1)
        assignment = result.assignments[0]
        self.assertEqual(assignment.molecule_index, 0)
        self.assertEqual(assignment.conformer_index, 0)
        self.assertEqual(assignment.source.record_id, "water-record")
        self.assertEqual(assignment.atom_ids, ("O", "H1", "H2"))
        self.assertEqual(assignment.conformer_id, "model-a")
        self.assertEqual(assignment.values.dtype, np.dtype(np.float64))
        self.assertTrue(assignment.values.flags.c_contiguous)
        self.assertFalse(assignment.values.flags.writeable)
        self.assertIsNotNone(assignment.values.base)
        self.assertTrue(np.isclose(assignment.values.sum(), 0.0))
        self.assertIs(result.requested, self.full_eem)
        effective = result.effective
        if effective is None:
            self.fail("successful calculation must report effective provenance")
        self.assertEqual(effective.method_id, "eem")
        self.assertEqual(effective.parameter_set_id, selected_candidate.parameter_set_id)
        self.assertIs(effective.execution_policy.mode, chargefw.ExecutionMode.FULL)
        self.assertGreaterEqual(result.timings.applicability_seconds, 0.0)
        self.assertGreaterEqual(result.timings.computation_seconds, 0.0)
        with self.assertRaises(TypeError):
            cast(Any, self.full_eem.method_options)["eem"] = {"unexpected": True}
        with self.assertRaises(RuntimeError):
            assessment.calculate()

    def test_assignment_cardinality_and_mapping(self) -> None:
        geometry_independent = self.calculator.calculate(
            chargefw.Molecule([8, 1, 1]),
            chargefw.CalculationOptions(
                method="formal", execution=chargefw.ExecutionSelectionKind.FULL
            ),
        )
        self.assertIs(geometry_independent.status, chargefw.ExecutionStatus.SUCCESS)
        self.assertEqual(len(geometry_independent.assignments), 1)
        self.assertIsNone(geometry_independent.assignments[0].conformer_index)

        collection = chargefw.MoleculeCollection([water(2)])
        options = chargefw.CalculationOptions(
            method="qeq", execution=chargefw.ExecutionSelectionKind.FULL
        )
        multi_result = self.calculator.calculate(collection, options)
        self.assertIs(multi_result.status, chargefw.ExecutionStatus.SUCCESS)
        self.assertEqual(
            [item.conformer_index for item in multi_result.assignments], [0, 1]
        )
        self.assertEqual(
            [item.conformer_id for item in multi_result.assignments],
            ["model-a", "model-b"],
        )
        repeated_result = self.calculator.calculate(collection, options)
        self.assertEqual(
            [item.values.tolist() for item in repeated_result.assignments],
            [item.values.tolist() for item in multi_result.assignments],
        )

    def test_reduced_execution_policies(self) -> None:
        reduced = self.calculator.calculate(
            water(),
            chargefw.CalculationOptions(
                method="eem",
                execution=chargefw.ExecutionSelectionKind.CUTOFF,
                radius=8.0,
                charge_correction=chargefw.ChargeCorrectionPolicy.NONE,
            ),
        )
        self.assertIs(reduced.status, chargefw.ExecutionStatus.SUCCESS)
        reduced_effective = reduced.effective
        if reduced_effective is None:
            self.fail("successful calculation must report effective provenance")
        self.assertEqual(
            reduced_effective.execution_policy,
            chargefw.ExecutionPolicy(
                mode=chargefw.ExecutionMode.CUTOFF,
                radius=8.0,
                charge_correction=chargefw.ChargeCorrectionPolicy.NONE,
            ),
        )

        covered = self.calculator.calculate(
            water(),
            chargefw.CalculationOptions(
                method="eem", execution=chargefw.ExecutionSelectionKind.COVER, radius=8.0
            ),
        )
        self.assertIs(covered.status, chargefw.ExecutionStatus.SUCCESS)
        covered_effective = covered.effective
        if covered_effective is None:
            self.fail("successful calculation must report effective provenance")
        self.assertIs(covered_effective.execution_policy.mode, chargefw.ExecutionMode.COVER)

    def test_no_plan_result_and_typed_exception(self) -> None:
        result = self.calculator.assess(
            chargefw.Molecule([8, 1, 1], bonds=[[0, 1, 1], [0, 2, 1]]),
            chargefw.CalculationOptions(
                method="qeq", execution=chargefw.ExecutionSelectionKind.FULL
            ),
        ).calculate()
        self.assertIs(result.status, chargefw.ExecutionStatus.NO_EXECUTABLE_PLAN)
        self.assertEqual(result.assignments, ())
        self.assertTrue(result.applicability.rejected)
        self.assertIsInstance(
            result.applicability.rejected[0].issues[0].kind,
            chargefw.PrerequisiteIssueKind,
        )
        with self.assertRaises(chargefw.NoExecutablePlanError) as context:
            result.raise_for_status()
        self.assertIs(context.exception.result, result)

    def test_invalid_options_are_rejected_early(self) -> None:
        invalid_options = (
            (
                ValueError,
                dict(
                    execution=chargefw.ExecutionSelectionKind.FULL,
                    charge_correction=chargefw.ChargeCorrectionPolicy.UNIFORM,
                ),
            ),
            (TypeError, dict(execution="full")),
            (TypeError, dict(charge_correction="none")),
            (TypeError, dict(max_threads=True)),
            (TypeError, dict(cutoff_atom_threshold=False)),
            (ValueError, dict(cover_atom_threshold=np.iinfo(np.uintp).max + 1)),
            (ValueError, dict(max_threads=np.iinfo(np.int32).max + 1)),
        )
        for error_type, options in invalid_options:
            with self.subTest(error_type=error_type, options=options):
                with self.assertRaises(error_type):
                    chargefw.CalculationOptions(**options)

    def test_catalogs_and_descriptors_are_immutable_values(self) -> None:
        other_calculator = chargefw.Calculator()
        self.assertIs(other_calculator._catalog, self.calculator._catalog)
        self.assertIs(other_calculator.parameter_sets, self.calculator.parameter_sets)
        self.assertIs(chargefw.method_descriptors(), chargefw.method_descriptors())

        methods = self.calculator.methods
        eem = next(method for method in methods if method.id == "eem")
        self.assertTrue(eem.requires_coordinates)
        self.assertTrue(eem.supports_cutoff)
        self.assertTrue(eem.supports_cover)
        self.assertIsInstance(eem.options, tuple)
        self.assertTrue(
            all(option.type is chargefw.MethodOptionType.INTEGER for option in eem.options)
        )
        self.assertEqual(chargefw.method_descriptors(), methods)

        parameter_directory = Path(chargefw.__file__).parent / "_data" / "parameters"
        eem_path = next(parameter_directory.glob("EEM_Ouy2009.json"))
        parameter_set = chargefw.load_parameter_set(eem_path)
        self.assertEqual(parameter_set.descriptor.method_id, "eem")
        self.assertEqual(parameter_set.id, parameter_set.descriptor.id)
        self.assertEqual(
            repr(parameter_set),
            f"ParameterSet(id={parameter_set.id!r}, method_id='eem')",
        )
        with self.assertRaises(TypeError):
            chargefw.ParameterSet()
        with self.assertRaises(AttributeError):
            parameter_set._native = None
        loaded_sets = chargefw.load_parameter_sets(parameter_directory)
        self.assertIn(parameter_set.id, {value.id for value in loaded_sets})

        explicit_calculator = chargefw.Calculator([parameter_set])
        self.assertEqual(explicit_calculator.parameter_sets, (parameter_set.descriptor,))
        self.assertEqual(repr(explicit_calculator), "Calculator(parameter_sets=1)")
        result = explicit_calculator.calculate(water(), self.full_eem)
        self.assertIs(result.status, chargefw.ExecutionStatus.SUCCESS)
        effective = result.effective
        if effective is None:
            self.fail("successful calculation must report effective provenance")
        self.assertEqual(effective.parameter_set_id, parameter_set.id)
        with self.assertRaises(ValueError):
            chargefw.Calculator([])


if __name__ == "__main__":
    unittest.main()
