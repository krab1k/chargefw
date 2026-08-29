"""Focused calculation facade and result-model checks."""

import gc
import unittest
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
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


def calculate_formal(molecules: Any) -> chargefw.CalculationResult:
    return chargefw.calculate(molecules, method="formal", execution="full", threads=1)


def assess_formal(molecules: Any) -> chargefw.Assessment:
    return chargefw.assess(molecules, method="formal", execution="full", threads=1)


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
    def test_default_calculation_selects_a_supported_plan(self) -> None:
        result = chargefw.calculate(water())
        self.assertEqual(result.status, "success")
        self.assertIsNotNone(result.selected)
        self.assertIsNotNone(result.effective)

    def test_assessment_and_full_calculation(self) -> None:
        assessment = chargefw.assess(water(), method="eem", execution="full")
        self.assertTrue(assessment.executable)
        self.assertEqual(
            assessment.execution,
            chargefw.ExecutionPolicy(
                mode="full",
                radius=None,
                charge_correction="none",
            ),
        )
        selected_candidate = assessment.selected
        if selected_candidate is None:
            self.fail("executable assessment must select a method")
        self.assertEqual(selected_candidate.method.id, "eem")
        self.assertEqual(selected_candidate.executions[0].mode, "full")

        result = assessment.calculate()
        self.assertEqual(result.status, "success")
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
        self.assertIs(np.asarray(assignment), assignment.values)
        self.assertEqual(np.asarray(assignment, dtype=np.float32).dtype, np.dtype(np.float32))
        self.assertEqual(result.requested.method, "eem")
        self.assertEqual(result.requested.execution, "full")
        effective = result.effective
        if effective is None:
            self.fail("successful calculation must report effective provenance")
        self.assertEqual(effective.method.id, "eem")
        self.assertEqual(effective.parameter_set, selected_candidate.parameter_set)
        self.assertEqual(effective.execution.mode, "full")
        self.assertGreaterEqual(result.timings.applicability_seconds, 0.0)
        self.assertGreaterEqual(result.timings.computation_seconds, 0.0)
        with self.assertRaises(TypeError):
            cast(Any, result.requested.options_by_method)["eem"] = {"unexpected": True}
        with self.assertRaises(RuntimeError):
            assessment.calculate()

    def test_assessment_context_manager_releases_prepared_state(self) -> None:
        with chargefw.assess(water(), method="eem", execution="full") as assessment:
            self.assertTrue(assessment.executable)
        with self.assertRaisesRegex(RuntimeError, "closed"):
            assessment.calculate()

    def test_assignment_cardinality_and_mapping(self) -> None:
        geometry_independent = chargefw.calculate(
            chargefw.Molecule([8, 1, 1]),
            method="formal",
            execution="full",
        )
        self.assertEqual(geometry_independent.status, "success")
        self.assertEqual(len(geometry_independent.assignments), 1)
        self.assertIsNone(geometry_independent.assignments[0].conformer_index)

        collection = chargefw.MoleculeCollection([water(2)])
        multi_result = chargefw.calculate(collection, method="qeq", execution="full")
        self.assertEqual(multi_result.status, "success")
        self.assertEqual([item.conformer_index for item in multi_result.assignments], [0, 1])
        self.assertEqual(
            [item.conformer_id for item in multi_result.assignments],
            ["model-a", "model-b"],
        )
        repeated_result = chargefw.calculate(collection, method="qeq", execution="full")
        self.assertEqual(
            [item.values.tolist() for item in repeated_result.assignments],
            [item.values.tolist() for item in multi_result.assignments],
        )

        mapped_collection = chargefw.MoleculeCollection(
            [
                chargefw.Molecule([1], source_name="first", record_id="record-a", atom_ids=["A"]),
                chargefw.Molecule(
                    [8, 1],
                    source_name="second",
                    record_id="record-b",
                    atom_ids=["B", "C"],
                ),
            ]
        )
        mapped_result = calculate_formal(mapped_collection)
        self.assertEqual(mapped_result.status, "success")
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
            with (
                self.subTest(error_type=error_type, overrides=overrides),
                self.assertRaises(error_type),
            ):
                chargefw.ChargeAssignment(**(defaults | overrides))

    def test_functional_api_supports_concurrent_calculations(self) -> None:
        collection = chargefw.MoleculeCollection([chargefw.Molecule([8, 1, 1])])
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = [executor.submit(calculate_formal, collection) for _ in range(4)]
        results = [future.result() for future in futures]
        self.assertTrue(all(result.status == "success" for result in results))
        self.assertEqual(
            [result.assignments[0].values.tolist() for result in results],
            [[0.0, 0.0, 0.0]] * 4,
        )

    def test_assessment_releases_the_gil(self) -> None:
        molecule = chargefw.Molecule([1] * _GIL_TEST_ATOM_COUNT)
        assessment, progressed, operation_seconds = run_while_python_thread_progresses(
            lambda: assess_formal(molecule)
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
        assessment = assess_formal(chargefw.Molecule([1] * _GIL_TEST_ATOM_COUNT))
        result, progressed, operation_seconds = run_while_python_thread_progresses(
            assessment.calculate
        )
        self.assertEqual(result.status, "success")
        self.assertGreater(operation_seconds, _GIL_OBSERVATION_DELAY_SECONDS)
        self.assertTrue(progressed)

    def test_reduced_execution_policies(self) -> None:
        reduced = chargefw.calculate(
            water(),
            method="eem",
            execution="cutoff",
            radius=8.0,
            charge_correction="none",
        )
        self.assertEqual(reduced.status, "success")
        reduced_effective = reduced.effective
        if reduced_effective is None:
            self.fail("successful calculation must report effective provenance")
        self.assertEqual(
            reduced_effective.execution,
            chargefw.ExecutionPolicy(
                mode="cutoff",
                radius=8.0,
                charge_correction="none",
            ),
        )

        covered = chargefw.calculate(
            water(),
            method="eem",
            execution="cover",
            radius=8.0,
        )
        self.assertEqual(covered.status, "success")
        covered_effective = covered.effective
        if covered_effective is None:
            self.fail("successful calculation must report effective provenance")
        self.assertEqual(covered_effective.execution.mode, "cover")

    def test_automatic_thresholds_and_explicit_full_warnings(self) -> None:
        automatic = chargefw.assess(
            water(),
            method="eem",
            cutoff_threshold=1,
            cover_threshold=100,
        )
        self.assertTrue(automatic.executable)
        if automatic.execution is None:
            self.fail("executable assessment must report an execution policy")
        self.assertEqual(automatic.execution.mode, "cutoff")
        self.assertEqual(automatic.execution.radius, 12.0)

        explicit_full = chargefw.assess(
            water(),
            method="eem",
            execution="full",
            cutoff_threshold=1,
            cover_threshold=100,
        )
        self.assertTrue(explicit_full.executable)
        self.assertTrue(explicit_full.warnings)
        self.assertEqual(explicit_full.warnings[0].kind, "resource_threshold_exceeded")

    def test_no_plan_result_and_typed_exception(self) -> None:
        assessment = chargefw.assess(
            chargefw.Molecule([8, 1, 1], bonds=[[0, 1, 1], [0, 2, 1]]),
            method="qeq",
            execution="full",
        )
        self.assertFalse(assessment.executable)
        with self.assertRaises(chargefw.NoExecutablePlanError) as context:
            assessment.calculate()
        result = context.exception.result
        self.assertEqual(result.status, "no_executable_plan")
        self.assertEqual(result.assignments, ())
        self.assertTrue(result.rejected)
        self.assertEqual(result.rejected[0].issues[0].kind, "missing_feature")

    def test_invalid_selection_requests_raise_value_error(self) -> None:
        invalid_requests: tuple[dict[str, Any], ...] = (
            {"method": "not-a-method"},
            {"method": "eem", "parameter_set": "not-a-parameter-set"},
        )
        for arguments in invalid_requests:
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                chargefw.assess(water(), **arguments)

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
            assessment = assess_formal(molecule)
            report = assessment.report
            return assessment.calculate(), report

        result, report = calculate_owned_result()
        gc.collect()
        self.assertTrue(report.executable)
        self.assertEqual(result.status, "success")
        self.assertEqual(result.assignments[0].source.record_id, "owned-record")
        self.assertEqual(result.assignments[0].atom_ids, ("O", "H1", "H2"))
        np.testing.assert_array_equal(result.assignments[0].values, [0.0, 0.0, 0.0])

    def test_invalid_options_are_rejected_early(self) -> None:
        invalid_options = (
            (ValueError, {"execution": "fast"}),
            (ValueError, {"charge_correction": "invalid"}),
            (ValueError, {"parameter_matching": "guess"}),
            (TypeError, {"threads": True}),
            (TypeError, {"cutoff_threshold": False}),
            (ValueError, {"cover_threshold": np.iinfo(np.uintp).max + 1}),
            (ValueError, {"threads": np.iinfo(np.int32).max + 1}),
        )
        for error_type, options in invalid_options:
            with (
                self.subTest(error_type=error_type, options=options),
                self.assertRaises(error_type),
            ):
                chargefw.assess(water(), **options)

        with self.assertRaisesRegex(ValueError, "requires an explicit method"):
            chargefw.assess(water(), options={"iters": 2})
        with self.assertRaisesRegex(ValueError, "cannot be used together"):
            chargefw.assess(
                water(),
                method="peoe",
                options={"iters": 2},
                options_by_method={"peoe": {"iters": 2}},
            )

    def test_method_option_overrides_are_validated_and_reported(self) -> None:
        molecule = chargefw.Molecule([8, 1, 1], bonds=[[0, 1, 1], [0, 2, 1]])
        result = chargefw.calculate(
            molecule,
            method="peoe",
            options={"iters": 2},
            execution="full",
        )
        self.assertEqual(result.status, "success")
        if result.effective is None:
            self.fail("successful calculation must report effective provenance")
        self.assertEqual(dict(result.effective.options), {"iters": 2})

        for options in (
            {"iters": 0},
            {"unknown": 1},
        ):
            with self.subTest(options=options), self.assertRaises(ValueError):
                chargefw.assess(
                    molecule,
                    method="peoe",
                    options=options,
                    execution="full",
                )
        with self.assertRaises(ValueError):
            chargefw.assess(
                molecule,
                method="peoe",
                options_by_method={"qeq": {"overlap_term": "Ohno"}},
                execution="full",
            )

    def test_catalogs_and_descriptors_are_immutable_values(self) -> None:
        self.assertFalse(hasattr(chargefw, "Calculator"))
        self.assertFalse(hasattr(chargefw, "load_parameter_set"))
        self.assertFalse(hasattr(chargefw, "load_parameter_sets"))
        self.assertFalse(hasattr(chargefw, "method_descriptors"))

        methods = chargefw.methods
        self.assertIsInstance(methods, chargefw.MethodCatalog)
        self.assertIn("eem", methods)
        self.assertEqual(methods.get("eem"), methods["eem"])
        self.assertIsNone(methods.get("not-a-method"))
        self.assertEqual(tuple(methods), tuple(method.id for method in methods.values()))
        self.assertEqual(methods["eem"].id, "eem")
        with self.assertRaisesRegex(KeyError, "unknown method ID"):
            methods["not-a-method"]
        with self.assertRaises(AttributeError):
            cast(Any, methods)._values = ()

        eem = methods["eem"]
        self.assertTrue(eem.requires_coordinates)
        self.assertTrue(eem.supports_cutoff)
        self.assertTrue(eem.supports_cover)
        self.assertEqual(
            tuple(eem.parameter_sets),
            tuple(chargefw.parameter_sets.for_method("eem")),
        )
        peoe = methods["peoe"]
        self.assertEqual(len(peoe.options), 1)
        self.assertIsInstance(peoe.options, chargefw.MethodOptionCatalog)
        self.assertEqual(peoe.options["iters"].id, "iters")
        self.assertEqual(peoe.options["iters"].type, "integer")
        self.assertEqual(peoe.options["iters"].default, 6)
        self.assertEqual(peoe.options["iters"].minimum, 1)
        qeq = methods["qeq"]
        overlap = qeq.options["overlap_term"]
        self.assertEqual(overlap.type, "string")
        self.assertIn("Ohno", overlap.choices)

        self.assertIsInstance(chargefw.parameter_sets, chargefw.ParameterSetCatalog)
        parameter_set = next(iter(eem.parameter_sets.values()))
        self.assertEqual(chargefw.parameter_sets[parameter_set.id], parameter_set)
        result = chargefw.calculate(
            water(),
            method=eem,
            parameter_set=parameter_set,
            execution="full",
        )
        self.assertEqual(result.status, "success")
        effective = result.effective
        if effective is None:
            self.fail("successful calculation must report effective provenance")
        self.assertEqual(effective.parameter_set, parameter_set)


if __name__ == "__main__":
    unittest.main()
