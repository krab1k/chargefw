"""Focused calculation facade and result-model checks."""

import gc
import json
import unittest
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event, Thread
from time import perf_counter, sleep
from typing import Any, TypeVar, cast

import chargefw
import numpy as np
from chargefw._chargefw import core as _native_core

T = TypeVar("T")
_GIL_OBSERVATION_DELAY_SECONDS = 0.001
_GIL_TEST_ATOM_COUNT = 500_000


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


def formal_options() -> chargefw.CalculationOptions:
    return chargefw.CalculationOptions(
        method="formal",
        execution=chargefw.ExecutionSelectionKind.FULL,
        max_threads=1,
    )


def run_while_python_thread_progresses(operation: Callable[[], T]) -> tuple[T, bool, float]:
    started = Event()
    begin = Event()
    observing = Event()
    progressed = Event()

    def observe() -> None:
        started.set()
        begin.wait()
        observing.set()
        sleep(_GIL_OBSERVATION_DELAY_SECONDS)
        progressed.set()

    worker = Thread(target=observe)
    worker.start()
    if not started.wait(timeout=1.0):
        worker.join(timeout=1.0)
        raise AssertionError("worker thread did not start")
    try:
        begin.set()
        if not observing.wait(timeout=1.0):
            raise AssertionError("worker thread did not begin observing")
        start = perf_counter()
        result = operation()
        operation_seconds = perf_counter() - start
        progressed_during_operation = progressed.is_set()
    finally:
        worker.join(timeout=1.0)
    if worker.is_alive():
        raise AssertionError("worker thread did not stop")
    return result, progressed_during_operation, operation_seconds


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
        with self.assertRaises(ValueError):
            assignment.values.setflags(write=True)
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

        mapped_collection = chargefw.MoleculeCollection(
            [
                chargefw.Molecule(
                    [1], source_name="first", record_id="record-a", atom_ids=["A"]
                ),
                chargefw.Molecule(
                    [8, 1],
                    source_name="second",
                    record_id="record-b",
                    atom_ids=["B", "C"],
                ),
            ]
        )
        mapped_result = self.calculator.calculate(mapped_collection, formal_options())
        self.assertIs(mapped_result.status, chargefw.ExecutionStatus.SUCCESS)
        self.assertEqual(
            [assignment.molecule_index for assignment in mapped_result.assignments], [0, 1]
        )
        self.assertEqual(
            [assignment.source.record_id for assignment in mapped_result.assignments],
            ["record-a", "record-b"],
        )
        self.assertEqual(
            [assignment.atom_ids for assignment in mapped_result.assignments],
            [("A",), ("B", "C")],
        )

    def test_public_charge_assignment_is_validated_and_comparable(self) -> None:
        source = chargefw.SourceIdentity("fixture", 1, "record")
        assignment = chargefw.ChargeAssignment(
            values=np.array([0.25, -0.25]),
            molecule_index=cast(Any, np.int64(2)),
            conformer_index=0,
            source=source,
            atom_ids=("A", "B"),
            conformer_id="model",
        )
        same = chargefw.ChargeAssignment(
            values=np.array([0.25, -0.25]),
            molecule_index=2,
            conformer_index=0,
            source=source,
            atom_ids=("A", "B"),
            conformer_id="model",
        )
        different = chargefw.ChargeAssignment(
            values=np.array([0.5, -0.5]),
            molecule_index=2,
            conformer_index=0,
            source=source,
            atom_ids=("A", "B"),
            conformer_id="model",
        )
        self.assertEqual(assignment, same)
        self.assertNotEqual(assignment, different)
        self.assertFalse(assignment.values.flags.writeable)
        with self.assertRaises(ValueError):
            assignment.values.setflags(write=True)

        invalid_assignments = (
            (ValueError, {"values": [[0.0]], "atom_ids": ("A",)}),
            (ValueError, {"values": [np.nan], "atom_ids": ("A",)}),
            (ValueError, {"values": [0.0], "atom_ids": ()}),
            (
                ValueError,
                {"values": [0.0], "atom_ids": ("A",), "molecule_index": -1},
            ),
            (
                ValueError,
                {"values": [0.0], "atom_ids": ("A",), "conformer_id": "model"},
            ),
            (
                TypeError,
                {"values": [0.0], "atom_ids": ("A",), "source": "fixture"},
            ),
        )
        defaults: dict[str, Any] = {
            "molecule_index": 0,
            "conformer_index": None,
            "source": source,
            "conformer_id": None,
        }
        for error_type, overrides in invalid_assignments:
            with self.subTest(error_type=error_type, overrides=overrides), self.assertRaises(
                error_type
            ):
                chargefw.ChargeAssignment(**(defaults | overrides))

    def test_shared_calculator_supports_concurrent_calculations(self) -> None:
        collection = chargefw.MoleculeCollection([chargefw.Molecule([8, 1, 1])])
        options = formal_options()
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = [
                executor.submit(self.calculator.calculate, collection, options) for _ in range(4)
            ]
        results = [future.result() for future in futures]
        self.assertTrue(
            all(result.status is chargefw.ExecutionStatus.SUCCESS for result in results)
        )
        self.assertEqual(
            [result.assignments[0].values.tolist() for result in results],
            [[0.0, 0.0, 0.0]] * 4,
        )

    def test_assessment_releases_the_gil(self) -> None:
        calculator = chargefw.Calculator()
        molecule = chargefw.Molecule([1] * _GIL_TEST_ATOM_COUNT)
        assessment, progressed, operation_seconds = run_while_python_thread_progresses(
            lambda: calculator.assess(molecule, formal_options())
        )
        self.assertTrue(assessment.executable)
        self.assertGreater(operation_seconds, _GIL_OBSERVATION_DELAY_SECONDS)
        self.assertTrue(progressed)

    def test_native_molecule_construction_releases_the_gil(self) -> None:
        atom_count = _GIL_TEST_ATOM_COUNT
        atomic_numbers = np.ones(atom_count, dtype=np.int64)
        formal_charges = np.zeros(atom_count, dtype=np.int64)
        bonds = np.empty((0, 3), dtype=np.int64)
        coordinates = np.empty((0, atom_count, 3), dtype=np.float64)
        atom_names = ("",) * atom_count
        molecule, progressed, operation_seconds = run_while_python_thread_progresses(
            lambda: _native_core._make_molecule(
                atomic_numbers,
                formal_charges,
                bonds,
                coordinates,
                atom_names,
                (),
                "gil-test",
            )
        )
        self.assertEqual(molecule.atom_count, atom_count)
        self.assertGreater(operation_seconds, _GIL_OBSERVATION_DELAY_SECONDS)
        self.assertTrue(progressed)

    def test_execution_releases_the_gil(self) -> None:
        calculator = chargefw.Calculator()
        assessment = calculator.assess(
            chargefw.Molecule([1] * _GIL_TEST_ATOM_COUNT), formal_options()
        )
        result, progressed, operation_seconds = run_while_python_thread_progresses(
            assessment.calculate
        )
        self.assertIs(result.status, chargefw.ExecutionStatus.SUCCESS)
        self.assertGreater(operation_seconds, _GIL_OBSERVATION_DELAY_SECONDS)
        self.assertTrue(progressed)

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

    def test_automatic_thresholds_and_explicit_full_warnings(self) -> None:
        automatic = self.calculator.assess(
            water(),
            chargefw.CalculationOptions(
                method="eem",
                cutoff_atom_threshold=1,
                cover_atom_threshold=100,
            ),
        )
        self.assertTrue(automatic.executable)
        if automatic.execution_policy is None:
            self.fail("executable assessment must report an execution policy")
        self.assertIs(automatic.execution_policy.mode, chargefw.ExecutionMode.CUTOFF)
        self.assertEqual(automatic.execution_policy.radius, 12.0)

        explicit_full = self.calculator.assess(
            water(),
            chargefw.CalculationOptions(
                method="eem",
                execution=chargefw.ExecutionSelectionKind.FULL,
                cutoff_atom_threshold=1,
                cover_atom_threshold=100,
            ),
        )
        self.assertTrue(explicit_full.executable)
        self.assertTrue(explicit_full.execution_issues)
        self.assertIs(
            explicit_full.execution_issues[0].kind,
            chargefw.ExecutionIssueKind.RESOURCE_THRESHOLD_EXCEEDED,
        )

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

    def test_numerical_failure_result_and_typed_exception(self) -> None:
        parameter_document = {
            "metadata": {"name": "Singular EEM", "method": "eem"},
            "common": {"names": ["kappa"], "values": [1.0]},
            "atom": {
                "names": ["A", "B"],
                "data": [
                    {"key": ["H", "plain", "*"], "value": [1.0, 1.0]},
                    {"key": ["O", "plain", "*"], "value": [2.0, 1.0]},
                ],
            },
        }
        with TemporaryDirectory() as directory:
            path = Path(directory) / "singular-eem.json"
            path.write_text(json.dumps(parameter_document), encoding="utf-8")
            calculator = chargefw.Calculator([chargefw.load_parameter_set(path)])
            result = calculator.calculate(
                chargefw.Molecule(
                    [1, 8],
                    coordinates=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                    name="singular-eem",
                ),
                chargefw.CalculationOptions(
                    method="eem",
                    parameter_set="singular-eem",
                    execution=chargefw.ExecutionSelectionKind.FULL,
                ),
            )
        self.assertIs(result.status, chargefw.ExecutionStatus.NUMERICAL_FAILURE)
        self.assertEqual(result.assignments, ())
        self.assertIsNotNone(result.effective)
        self.assertIn(
            "method 'eem', molecule 1 ('singular-eem'), conformer 1",
            result.failure_message or "",
        )
        with self.assertRaises(chargefw.NumericalFailureError) as context:
            result.raise_for_status()
        self.assertIs(context.exception.result, result)

    def test_invalid_selection_requests_raise_value_error(self) -> None:
        for options in (
            chargefw.CalculationOptions(method="not-a-method"),
            chargefw.CalculationOptions(method="eem", parameter_set="not-a-parameter-set"),
        ):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.calculator.assess(water(), options)

    def test_results_outlive_calculation_inputs(self) -> None:
        def calculate_owned_result() -> tuple[
            chargefw.CalculationResult, chargefw.AssessmentReport
        ]:
            molecule = chargefw.Molecule(
                [8, 1, 1],
                source_name="owned",
                record_id="owned-record",
                atom_ids=["O", "H1", "H2"],
            )
            calculator = chargefw.Calculator()
            assessment = calculator.assess(molecule, formal_options())
            report = assessment.report
            return assessment.calculate(), report

        result, report = calculate_owned_result()
        gc.collect()
        self.assertTrue(report.executable)
        self.assertIs(result.status, chargefw.ExecutionStatus.SUCCESS)
        self.assertEqual(result.assignments[0].source.record_id, "owned-record")
        self.assertEqual(result.assignments[0].atom_ids, ("O", "H1", "H2"))
        np.testing.assert_array_equal(result.assignments[0].values, [0.0, 0.0, 0.0])

    def test_invalid_options_are_rejected_early(self) -> None:
        invalid_options = (
            (
                ValueError,
                {
                    "execution": chargefw.ExecutionSelectionKind.FULL,
                    "charge_correction": chargefw.ChargeCorrectionPolicy.UNIFORM,
                },
            ),
            (TypeError, {"execution": "full"}),
            (TypeError, {"charge_correction": "none"}),
            (TypeError, {"max_threads": True}),
            (TypeError, {"cutoff_atom_threshold": False}),
            (ValueError, {"cover_atom_threshold": np.iinfo(np.uintp).max + 1}),
            (ValueError, {"max_threads": np.iinfo(np.int32).max + 1}),
        )
        for error_type, options in invalid_options:
            with self.subTest(error_type=error_type, options=options), self.assertRaises(
                error_type
            ):
                chargefw.CalculationOptions(**options)

    def test_method_option_overrides_are_validated_and_reported(self) -> None:
        molecule = chargefw.Molecule([8, 1, 1], bonds=[[0, 1, 1], [0, 2, 1]])
        options = chargefw.CalculationOptions(
            method="peoe",
            method_options={"peoe": {"iters": 2}},
            execution=chargefw.ExecutionSelectionKind.FULL,
        )
        result = self.calculator.calculate(molecule, options)
        self.assertIs(result.status, chargefw.ExecutionStatus.SUCCESS)
        if result.effective is None:
            self.fail("successful calculation must report effective provenance")
        self.assertEqual(dict(result.effective.method_options), {"iters": 2})

        for method_options in (
            {"peoe": {"iters": 0}},
            {"peoe": {"unknown": 1}},
            {"qeq": {"overlap_term": "Ohno"}},
        ):
            with self.subTest(method_options=method_options), self.assertRaises(ValueError):
                self.calculator.assess(
                    molecule,
                    chargefw.CalculationOptions(
                        method="peoe",
                        method_options=method_options,
                        execution=chargefw.ExecutionSelectionKind.FULL,
                    ),
                )

    def test_catalogs_and_descriptors_are_immutable_values(self) -> None:
        other_calculator = chargefw.Calculator()
        self.assertIs(other_calculator._catalog, self.calculator._catalog)
        self.assertIs(other_calculator.parameter_sets, self.calculator.parameter_sets)
        self.assertFalse(hasattr(chargefw, "method_descriptors"))

        methods = self.calculator.methods
        self.assertIsInstance(methods, chargefw.MethodCatalog)
        self.assertIn("eem", methods)
        self.assertEqual(methods.get("eem"), methods["eem"])
        self.assertIsNone(methods.get("not-a-method"))
        self.assertEqual(methods.ids(), tuple(method.id for method in methods))
        self.assertEqual(methods[:2], tuple(methods)[:2])
        with self.assertRaisesRegex(KeyError, "unknown method ID"):
            methods["not-a-method"]
        with self.assertRaises(AttributeError):
            cast(Any, methods)._values = ()

        eem = methods["eem"]
        self.assertTrue(eem.requires_coordinates)
        self.assertTrue(eem.supports_cutoff)
        self.assertTrue(eem.supports_cover)
        self.assertEqual(
            eem.parameter_sets.ids(),
            self.calculator.parameter_sets.for_method("eem").ids(),
        )
        peoe = methods["peoe"]
        self.assertEqual(len(peoe.options), 1)
        self.assertIsInstance(peoe.options, chargefw.MethodOptionCatalog)
        self.assertEqual(peoe.options["iters"].id, "iters")
        self.assertIs(peoe.options["iters"].type, chargefw.MethodOptionType.INTEGER)
        self.assertEqual(peoe.options["iters"].default, 6)
        self.assertEqual(peoe.options["iters"].minimum, 1)
        qeq = methods["qeq"]
        overlap = qeq.options["overlap_term"]
        self.assertIs(overlap.type, chargefw.MethodOptionType.STRING)
        self.assertIn("Ohno", overlap.choices)

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
            cast(Any, parameter_set)._native = None
        loaded_sets = chargefw.load_parameter_sets(parameter_directory)
        self.assertIn(parameter_set.id, {value.id for value in loaded_sets})

        explicit_calculator = chargefw.Calculator([parameter_set])
        self.assertIsInstance(explicit_calculator.parameter_sets, chargefw.ParameterSetCatalog)
        self.assertEqual(explicit_calculator.parameter_sets[0], parameter_set.descriptor)
        self.assertEqual(explicit_calculator.parameter_sets[parameter_set.id], parameter_set.descriptor)
        self.assertEqual(
            explicit_calculator.methods["eem"].parameter_sets[parameter_set.id],
            parameter_set.descriptor,
        )
        self.assertEqual(repr(explicit_calculator), "Calculator(parameter_sets=1)")
        result = explicit_calculator.calculate(
            water(),
            chargefw.CalculationOptions(
                method=explicit_calculator.methods["eem"],
                parameter_set=explicit_calculator.parameter_sets[parameter_set.id],
                execution=chargefw.ExecutionSelectionKind.FULL,
            ),
        )
        self.assertIs(result.status, chargefw.ExecutionStatus.SUCCESS)
        effective = result.effective
        if effective is None:
            self.fail("successful calculation must report effective provenance")
        self.assertEqual(effective.parameter_set_id, parameter_set.id)
        with self.assertRaises(ValueError):
            chargefw.Calculator([])
        with self.assertRaises(ValueError):
            chargefw.Calculator([parameter_set, parameter_set])


if __name__ == "__main__":
    unittest.main()
