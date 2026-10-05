"""Focused calculation facade and result-model checks."""

import gc
import unittest
from collections.abc import Callable, Mapping
from concurrent.futures import ThreadPoolExecutor
from threading import Event, Lock, Thread
from time import perf_counter, sleep
from typing import Any, TypeVar, cast
from unittest.mock import patch

import chargefw
import numpy as np
from chargefw._chargefw import calculation as _native_calculation
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
    )


FIXED_IONS_MMCIF = """data_fixed_groups
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1 O O . HOH W . ? 0 0 0 1 20 0 5 HOH W O 1
HETATM 2 H H1 . HOH W . ? 0.96 0 0 1 20 0 5 HOH W H1 1
HETATM 3 H H2 . HOH W . ? -0.24 0.93 0 1 20 0 5 HOH W H2 1
HETATM 4 Mg MG . MG I . ? 0 0 3 1 20 0 9 MG I MG 1
#
"""


def fixed_groups_molecules() -> chargefw.MoleculeCollection:
    return chargefw.io.parse(FIXED_IONS_MMCIF, format="mmcif", bonds="templates")


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
    def test_ion_selections_are_immutable_native_catalog_tuples(self) -> None:
        common: tuple[str, ...] = chargefw.COMMON_IONS
        all_ions: tuple[str, ...] = chargefw.ALL_IONS
        self.assertEqual(common, ("NA", "K", "MG", "CA", "CL", "ZN", "FE", "FE2"))
        self.assertIsInstance(all_ions, tuple)
        self.assertLess(set(common), set(all_ions))
        self.assertEqual(len(all_ions), len(set(all_ions)))
        self.assertEqual(len(all_ions), 80)
        self.assertEqual(all_ions, tuple(sorted(all_ions)))
        self.assertEqual(common, _native_calculation._fixed_ion_names(True))
        self.assertEqual(all_ions, _native_calculation._fixed_ion_names(False))
        for name in ("COMMON_IONS", "ALL_IONS"):
            self.assertIn(name, chargefw.__all__)
        for selection in (common, all_ions):
            with self.assertRaises(TypeError):
                cast(Any, selection)[0] = "MG"
            with self.assertRaises(AttributeError):
                cast(Any, selection).append("MG")
            absent = chargefw.calculate(
                chargefw.io.parse(FIXED_IONS_MMCIF.split("HETATM 4")[0] + "#\n", format="mmcif"),
                method="formal",
                fixed_ions=selection,
            )
            self.assertEqual(absent.requested.fixed_ions, selection)
            if absent.plan is None:
                self.fail("absent ion selection must retain the formal plan")
            self.assertIsNone(absent.plan.fixed_ions)

    def test_all_ions_selects_noncommon_template_in_assessment_and_calculation(self) -> None:
        # CU1's atom name is CU and its template charge differs from the imported zero.
        text = FIXED_IONS_MMCIF.replace(
            "HETATM 4 Mg MG . MG I . ? 0 0 3 1 20 0 9 MG I MG 1",
            "HETATM 4 Cu CU . CU1 I . ? 0 0 3 1 20 0 9 CU1 I CU 1",
        )
        molecules = chargefw.io.parse(text, format="mmcif", bonds="templates")
        common = chargefw.calculate(molecules, method="formal", fixed_ions=chargefw.COMMON_IONS)
        if common.plan is None:
            self.fail("unselected noncommon ion must leave formal calculation available")
        self.assertIsNone(common.plan.fixed_ions)
        np.testing.assert_array_equal(common.assignments[0].values, [0.0] * 4)
        assessment = chargefw.assess(
            molecules,
            method="sqeqp",
            parameter_set="SQEqp_Schindler2021_CCD_gen",
            execution="full",
            fixed_ions=chargefw.ALL_IONS,
        )
        if assessment.default_plan is None:
            self.fail("ALL_IONS must produce an SQEqp plan with the CU1 template")
        result = chargefw.calculate(molecules, assessment.default_plan)
        self.assertEqual(result.requested.fixed_ions, chargefw.ALL_IONS)
        self.assertEqual(result.assignments[0].values[3], 1.0)
        self.assertAlmostEqual(float(result.assignments[0].values[:3].sum()), 0.0, places=12)
        np.testing.assert_array_equal(molecules[0].formal_charges, [0] * 4)
        if result.plan is None or result.plan.fixed_ions is None:
            self.fail("ALL_IONS must retain effective source provenance")
        self.assertEqual(result.plan.fixed_ions.sources, (chargefw.FixedAtomCharge(0, 3, 1.0),))

    def test_observer_receives_owned_execution_progress(self) -> None:
        class RecordingObserver(chargefw.CalculationObserver):
            def __init__(self) -> None:
                self.events: list[chargefw.CalculationProgress] = []

            def on_progress(self, progress: chargefw.CalculationProgress) -> None:
                self.events.append(progress)

        observer = RecordingObserver()
        collection = chargefw.MoleculeCollection(
            [chargefw.Molecule([1]), chargefw.Molecule([8, 1])]
        )
        result = chargefw.calculate(
            collection,
            method="formal",
            execution="full",
            threads=1,
            observer=observer,
        )

        self.assertEqual(result.status, "success")
        self.assertEqual(
            [event.phase for event in observer.events],
            [
                "computation_started",
                "target_started",
                "target_finished",
                "target_started",
                "target_finished",
                "computation_finished",
            ],
        )
        started = observer.events[0]
        self.assertEqual(started.method_id, "formal")
        self.assertEqual(started.mode, "full")
        targets = [event for event in observer.events if event.phase == "target_started"]
        self.assertEqual([event.target_index for event in targets], [0, 1])
        self.assertEqual([event.target_count for event in targets], [2, 2])
        self.assertEqual([event.molecule_index for event in targets], [0, 1])
        self.assertEqual([event.conformer_index for event in targets], [None, None])
        self.assertTrue(all(event.elapsed_seconds >= 0.0 for event in observer.events))

    def test_observer_receives_reduced_fragment_progress(self) -> None:
        class RecordingObserver(chargefw.CalculationObserver):
            def __init__(self) -> None:
                self.events: list[chargefw.CalculationProgress] = []

            def on_progress(self, progress: chargefw.CalculationProgress) -> None:
                self.events.append(progress)

        observer = RecordingObserver()
        result = chargefw.calculate(
            water(),
            method="eem",
            execution="cutoff",
            radius=8.0,
            threads=1,
            observer=observer,
        )

        self.assertEqual(result.status, "success")
        fragments = [event for event in observer.events if event.phase == "fragment_progress"]
        self.assertTrue(fragments)
        self.assertEqual(fragments[-1].completed_fragment_count, fragments[-1].fragment_count)
        self.assertEqual(fragments[-1].target_index, 0)
        self.assertEqual(fragments[-1].target_count, 1)

    def test_fixed_ions_snapshot_resolve_and_reuse_plan(self) -> None:
        molecules = fixed_groups_molecules()
        self.assertEqual(molecules[0].bond_count, 2)
        names = ["MG"]
        assessment = chargefw.assess(
            molecules,
            method="sqeqp",
            parameter_set="SQEqp_Schindler2021_CCD_gen",
            execution="full",
            fixed_ions=names,
        )
        names.append("CA")
        plan = assessment.default_plan
        if plan is None:
            self.fail("assessment must produce an SQEqp plan with the selected MG source")
        del assessment
        gc.collect()
        result = chargefw.calculate(molecules, plan)
        direct = chargefw.calculate(
            molecules,
            method="sqeqp",
            parameter_set="SQEqp_Schindler2021_CCD_gen",
            execution="full",
            fixed_ions=chargefw.COMMON_IONS,
        )
        self.assertEqual(result.requested.fixed_ions, ("MG",))
        self.assertEqual(result.assignments[0].values[3], 2.0)
        np.testing.assert_allclose(result.assignments[0].values, direct.assignments[0].values)
        self.assertEqual(molecules[0].formal_charges.tolist(), [0, 0, 0, 0])
        if result.plan is None or result.plan.fixed_ions is None:
            self.fail("executed plan must retain effective source provenance")
        self.assertEqual(
            result.plan.fixed_ions.sources,
            (chargefw.FixedAtomCharge(0, 3, 2.0),),
        )
        self.assertEqual(
            result.plan.fixed_ions, chargefw.FixedIons([chargefw.FixedAtomCharge(0, 3, 2.0)])
        )
        self.assertFalse(hasattr(result.plan.fixed_ions, "charge_provenance"))
        self.assertFalse(hasattr(result.plan.fixed_ions, "charge_totals"))
        with self.assertRaises(AttributeError):
            cast(Any, result.plan.fixed_ions).sources = ()
        self.assertFalse(hasattr(chargefw, "FixedIonChargeTotals"))
        self.assertFalse(hasattr(chargefw, "FixedIonsProvenance"))
        with self.assertRaises(AttributeError):
            cast(Any, result.requested).fixed_ions = ("CA",)
        with self.assertRaisesRegex(TypeError, "selection arguments"):
            chargefw.calculate(molecules, plan, fixed_ions=[])

    def test_fixed_ions_validation_empty_and_method_support(self) -> None:
        molecules = fixed_groups_molecules()
        for invalid in ("MG", b"MG", ["MG", 12]):
            with self.subTest(invalid=invalid), self.assertRaises(TypeError):
                chargefw.assess(molecules, fixed_ions=cast(Any, invalid))
        with self.assertRaisesRegex(ValueError, "unknown fixed-ion component ID"):
            chargefw.assess(molecules, fixed_ions=["mg"])
        with self.assertRaisesRegex(ValueError, "unknown fixed-ion component ID"):
            chargefw.assess(chargefw.MoleculeCollection([]), fixed_ions=["bad"])

        absent = chargefw.calculate(molecules, method="formal", fixed_ions=["CA"])
        self.assertEqual(absent.requested.fixed_ions, ("CA",))
        if absent.plan is None:
            self.fail("formal calculation must produce a plan")
        self.assertIsNone(absent.plan.fixed_ions)
        duplicate = chargefw.calculate(molecules, method="eem", fixed_ions=["MG", "MG"])
        if duplicate.plan is None or duplicate.plan.fixed_ions is None:
            self.fail("duplicate names must resolve to one fixed source")
        self.assertEqual(len(duplicate.plan.fixed_ions.sources), 1)
        self.assertEqual(chargefw.RequestedCalculation().fixed_ions, ())

        plain = chargefw.calculate(water(), method="formal")
        empty = chargefw.calculate(water(), method="formal", fixed_ions=[])
        self.assertEqual(empty.requested.fixed_ions, ())
        np.testing.assert_array_equal(plain.assignments[0].values, empty.assignments[0].values)

        unsupported = chargefw.assess(
            molecules, method="qeq", parameter_set="QEq_original", fixed_ions=["MG"]
        )
        self.assertFalse(unsupported.plans)
        self.assertIn(
            "unsupported_fixed_ions",
            [issue.kind for rejection in unsupported.rejections for issue in rejection.issues],
        )
        with self.assertRaisesRegex(ValueError, "component metadata"):
            chargefw.assess(water(), fixed_ions=["MG"])

    def test_observer_receives_parallel_target_progress(self) -> None:
        class CountingObserver(chargefw.CalculationObserver):
            def __init__(self) -> None:
                self.lock = Lock()
                self.completed = 0

            def on_progress(self, progress: chargefw.CalculationProgress) -> None:
                if progress.phase == "target_finished":
                    with self.lock:
                        self.completed += 1

        molecule_count = 64
        observer = CountingObserver()
        result = chargefw.calculate(
            chargefw.MoleculeCollection(chargefw.Molecule([1]) for _ in range(molecule_count)),
            method="formal",
            execution="full",
            threads=2,
            observer=observer,
        )

        self.assertEqual(result.status, "success")
        self.assertEqual(observer.completed, molecule_count)

    def test_observer_can_cancel_a_reusable_plan(self) -> None:
        class CancellingObserver(chargefw.CalculationObserver):
            def __init__(self) -> None:
                self.cancel_requested = False
                self.phases: list[str] = []

            def on_progress(self, progress: chargefw.CalculationProgress) -> None:
                self.phases.append(progress.phase)
                if progress.phase == "target_started":
                    self.cancel_requested = True

            def cancelled(self) -> bool:
                return self.cancel_requested

        cases = (
            (chargefw.MoleculeCollection([water()]), None),
            (fixed_groups_molecules(), ["MG"]),
        )
        for molecule, groups in cases:
            with self.subTest(fixed_ions=groups):
                observer = CancellingObserver()
                options: dict[str, Any] = {
                    "method": "eem" if groups is None else "sqeqp",
                    "execution": "full",
                    "fixed_ions": groups,
                }
                if groups is not None:
                    options["parameter_set"] = "SQEqp_Schindler2021_CCD_gen"
                assessment = chargefw.assess(molecule, **options)
                plan = assessment.default_plan
                if plan is None:
                    self.fail("assessment must produce a default plan")
                with self.assertRaises(chargefw.CalculationCancelledError) as raised:
                    chargefw.calculate(molecule, plan, threads=1, observer=observer)

                cancelled = raised.exception.result
                self.assertEqual(cancelled.status, "cancelled")
                self.assertEqual(cancelled.assignments, ())
                self.assertEqual(observer.phases[0], "computation_started")
                self.assertIn("target_started", observer.phases)
                self.assertEqual(observer.phases[-1], "computation_finished")
                self.assertEqual(cancelled.requested.fixed_ions, tuple(groups or ()))
                if groups is not None:
                    self.assertIsNotNone(cancelled.plan)
                    assert cancelled.plan is not None
                    self.assertIsNotNone(cancelled.plan.fixed_ions)
                    assert cancelled.plan.fixed_ions is not None
                    self.assertEqual(
                        cancelled.plan.fixed_ions.sources,
                        (chargefw.FixedAtomCharge(0, 3, 2.0),),
                    )

                result = chargefw.calculate(molecule, plan, threads=1)
                self.assertEqual(result.status, "success")
                self.assertEqual(len(result.assignments), molecule[0].conformer_count)
                self.assertTrue(
                    all(
                        assignment.values.shape == (molecule[0].atom_count,)
                        for assignment in result.assignments
                    )
                )
                self.assertTrue(
                    all(np.all(np.isfinite(assignment.values)) for assignment in result.assignments)
                )

    def test_observer_must_use_public_base_class(self) -> None:
        with self.assertRaisesRegex(TypeError, "CalculationObserver"):
            chargefw.calculate(water(), observer=cast(Any, object()))

    def test_observer_callback_failure_does_not_change_calculation(self) -> None:
        class FailingObserver(chargefw.CalculationObserver):
            def on_progress(self, progress: chargefw.CalculationProgress) -> None:
                raise RuntimeError(f"failed during {progress.phase}")

        with patch("sys.unraisablehook") as unraisable:
            result = chargefw.calculate(
                water(),
                method="formal",
                execution="full",
                threads=1,
                observer=FailingObserver(),
            )

        self.assertEqual(result.status, "success")
        self.assertTrue(unraisable.called)
        # The recorded exception owns its traceback, which otherwise retains this calculation frame
        # until interpreter shutdown and produces a false-positive Nanobind leak report.
        unraisable.reset_mock()

    def test_default_calculation_selects_a_supported_plan(self) -> None:
        result = chargefw.calculate(water())
        self.assertEqual(result.status, "success")
        self.assertIsNotNone(result.plan)

    def test_assessment_and_full_calculation(self) -> None:
        molecule = water()
        assessment = chargefw.assess(molecule, method="eem", execution="full")
        self.assertTrue(assessment.plans)
        self.assertTrue(all(plan.policy.mode == "full" for plan in assessment.plans))
        plan = assessment.default_plan
        if plan is None:
            self.fail("assessment must produce a default plan")
        self.assertEqual(
            plan.policy,
            chargefw.ExecutionPolicy(
                mode="full",
                radius=None,
            ),
        )
        self.assertEqual(plan.method.id, "eem")

        result = chargefw.calculate(molecule, plan)
        self.assertEqual(result.status, "success")
        self.assertEqual(len(result.assignments), 1)
        assignment = result.assignments[0]
        self.assertEqual(assignment.molecule_index, 0)
        self.assertEqual(assignment.conformer_index, 0)
        self.assertEqual(assignment.source.record_id, "water-record")
        self.assertEqual(assignment.atom_ids, ("O", "H1", "H2"))
        self.assertEqual(assignment.values.dtype, np.dtype(np.float64))
        self.assertTrue(assignment.values.flags.c_contiguous)
        self.assertFalse(assignment.values.flags.writeable)
        self.assertEqual(assignment.values.ctypes.data % assignment.values.dtype.alignment, 0)
        with self.assertRaises(ValueError):
            assignment.values.setflags(write=True)
        self.assertTrue(np.isclose(assignment.values.sum(), 0.0))
        self.assertIs(np.asarray(assignment), assignment.values)
        self.assertEqual(np.asarray(assignment, dtype=np.float32).dtype, np.dtype(np.float32))
        self.assertEqual(result.requested.method, "eem")
        self.assertEqual(result.requested.execution, "full")
        executed = result.plan
        if executed is None:
            self.fail("successful calculation must report plan provenance")
        self.assertEqual(executed.method.id, "eem")
        self.assertEqual(executed.parameter_set, plan.parameter_set)
        self.assertEqual(executed.policy.mode, "full")
        self.assertGreaterEqual(result.timings.applicability_seconds, 0.0)
        self.assertGreaterEqual(result.timings.computation_seconds, 0.0)
        with self.assertRaises(TypeError):
            cast(Any, result.requested.options_by_method)["eem"] = {"unexpected": True}
        repeated = chargefw.calculate(molecule, plan)
        np.testing.assert_allclose(repeated.assignments[0].values, assignment.values)

        overridden_threads = chargefw.calculate(molecule, plan, threads=1)
        self.assertEqual(overridden_threads.requested.threads, 1)

    def test_plan_keeps_prepared_state_alive(self) -> None:
        molecule = water()
        plan = chargefw.assess(molecule, method="eem", execution="full").default_plan
        if plan is None:
            self.fail("assessment must produce a default plan")
        gc.collect()
        self.assertEqual(chargefw.calculate(molecule, plan).status, "success")

    def test_plans_are_filtered_reusable_and_target_bound(self) -> None:
        molecule = water()
        assessment = chargefw.assess(molecule, method="eqeq")
        self.assertEqual(
            {candidate_plan.policy.mode for candidate_plan in assessment.plans},
            {"full", "cutoff", "cover"},
        )
        for candidate_plan in assessment.plans:
            self.assertEqual(
                chargefw.calculate(molecule, candidate_plan, threads=1).status, "success"
            )

        cutoff = chargefw.assess(
            molecule,
            method="eqeq",
            execution="cutoff",
            radius=8.0,
        )
        self.assertTrue(cutoff.plans)
        self.assertTrue(all(plan.policy.mode == "cutoff" for plan in cutoff.plans))

        default_plan = assessment.default_plan
        if default_plan is None:
            self.fail("assessment must produce a default plan")
        with self.assertRaisesRegex(ValueError, "different molecule collection"):
            chargefw.calculate(water(), default_plan)
        with self.assertRaisesRegex(TypeError, "omit it"):
            chargefw.calculate(molecule, cast(Any, None))
        with self.assertRaisesRegex(TypeError, "cannot be combined"):
            chargefw.calculate(molecule, default_plan, method="eqeq")

        with ThreadPoolExecutor(max_workers=3) as executor:
            results = list(
                executor.map(
                    lambda _: chargefw.calculate(molecule, default_plan, threads=1),
                    range(3),
                )
            )
        self.assertTrue(all(result.status == "success" for result in results))

    def test_assignment_cardinality_and_mapping(self) -> None:
        geometry_independent = chargefw.calculate(
            chargefw.Molecule([8, 1, 1]),
            method="formal",
            execution="full",
        )
        self.assertEqual(geometry_independent.status, "success")
        self.assertEqual(len(geometry_independent.assignments), 1)
        self.assertIsNone(geometry_independent.assignments[0].conformer_index)
        self.assertIs(
            geometry_independent.assignment(molecule=0), geometry_independent.assignments[0]
        )
        with self.assertRaisesRegex(ValueError, "no conformer-specific"):
            geometry_independent.assignment(molecule=0, conformer=0)

        collection = chargefw.MoleculeCollection([water(2)])
        multi_result = chargefw.calculate(collection, method="qeq", execution="full")
        self.assertEqual(multi_result.status, "success")
        self.assertEqual([item.conformer_index for item in multi_result.assignments], [0, 1])
        self.assertEqual(multi_result.assignments_by_molecule, (multi_result.assignments,))
        self.assertIs(multi_result.assignment(molecule=0, conformer=0), multi_result.assignments[0])
        self.assertIs(multi_result.assignment(molecule=0, conformer=1), multi_result.assignments[1])
        with self.assertRaisesRegex(ValueError, "conformer is required"):
            multi_result.assignment(molecule=0)
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
        self.assertEqual(
            mapped_result.assignments_by_molecule,
            ((mapped_result.assignments[0],), (mapped_result.assignments[1],)),
        )

    def test_public_charge_assignment_is_validated_and_comparable(self) -> None:
        source = chargefw.SourceIdentity("fixture", 1, "record")
        assignment = chargefw.ChargeAssignment(
            values=np.array([0.25, -0.25]),
            molecule_index=cast(Any, np.int64(2)),
            conformer_index=0,
            source=source,
            atom_ids=("A", "B"),
        )
        same = chargefw.ChargeAssignment(
            values=np.array([0.25, -0.25]),
            molecule_index=2,
            conformer_index=0,
            source=source,
            atom_ids=("A", "B"),
        )
        different = chargefw.ChargeAssignment(
            values=np.array([0.5, -0.5]),
            molecule_index=2,
            conformer_index=0,
            source=source,
            atom_ids=("A", "B"),
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
                TypeError,
                {"values": [0.0], "atom_ids": ("A",), "source": "fixture"},
            ),
        )
        defaults: dict[str, Any] = {
            "molecule_index": 0,
            "conformer_index": None,
            "source": source,
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
        self.assertTrue(assessment.plans)
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
        self.assertIsInstance(molecule, _native_core._NativeMolecule)
        self.assertGreater(operation_seconds, _GIL_OBSERVATION_DELAY_SECONDS)
        self.assertTrue(progressed)

    def test_execution_releases_the_gil(self) -> None:
        molecule = chargefw.Molecule([1] * _GIL_TEST_ATOM_COUNT)
        assessment = assess_formal(molecule)
        plan = assessment.default_plan
        if plan is None:
            self.fail("assessment must produce a default plan")
        result, progressed, operation_seconds = run_while_python_thread_progresses(
            lambda: chargefw.calculate(molecule, plan)
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
        )
        self.assertEqual(reduced.status, "success")
        reduced_effective = reduced.plan
        if reduced_effective is None:
            self.fail("successful calculation must report plan provenance")
        self.assertEqual(
            reduced_effective.policy,
            chargefw.ExecutionPolicy(
                mode="cutoff",
                radius=8.0,
            ),
        )

        covered = chargefw.calculate(
            water(),
            method="eem",
            execution="cover",
            radius=8.0,
        )
        self.assertEqual(covered.status, "success")
        covered_effective = covered.plan
        if covered_effective is None:
            self.fail("successful calculation must report effective provenance")
        self.assertEqual(covered_effective.policy.mode, "cover")

    def test_automatic_thresholds_and_explicit_full_warnings(self) -> None:
        automatic = chargefw.assess(
            water(),
            method="eem",
            cutoff_threshold=1,
            cover_threshold=100,
        )
        if automatic.default_plan is None:
            self.fail("assessment must produce a default plan")
        self.assertEqual(automatic.default_plan.policy.mode, "cutoff")
        self.assertEqual(automatic.default_plan.policy.radius, 12.0)
        self.assertTrue(all(plan.policy.mode != "full" for plan in automatic.plans))
        self.assertTrue(
            any(
                rejection.policy is not None and rejection.policy.mode == "full"
                for rejection in automatic.rejections
            )
        )

        explicit_full = chargefw.assess(
            water(),
            method="eem",
            execution="full",
            cutoff_threshold=1,
            cover_threshold=100,
        )
        if explicit_full.default_plan is None:
            self.fail("explicit full assessment must produce a plan")
        self.assertTrue(all(plan.policy.mode == "full" for plan in explicit_full.plans))
        self.assertTrue(explicit_full.default_plan.warnings)
        self.assertEqual(explicit_full.default_plan.warnings[0].kind, "resource_threshold_exceeded")

    def test_no_plan_result_and_typed_exception(self) -> None:
        molecule = chargefw.Molecule([8, 1, 1], bonds=[[0, 1, 1], [0, 2, 1]])
        assessment = chargefw.assess(
            molecule,
            method="qeq",
            execution="full",
        )
        self.assertFalse(assessment.plans)
        with self.assertRaises(chargefw.NoExecutablePlanError) as context:
            chargefw.calculate(molecule, method="qeq", execution="full")
        result = context.exception.result
        self.assertEqual(result.status, "no_executable_plan")
        self.assertEqual(result.assignments, ())
        self.assertTrue(result.rejections)
        self.assertEqual(result.rejections[0].issues[0].kind, "missing_feature")

    def test_nonfinite_coordinates_are_assessed_per_method(self) -> None:
        molecule = chargefw.Molecule(
            [8, 1, 1],
            bonds=[[0, 1, 1], [0, 2, 1]],
            coordinates=[[0.0, 0.0, 0.0], [np.nan, 0.0, 0.0], [-0.24, 0.93, 0.0]],
        )
        self.assertTrue(np.isnan(molecule.coordinates[0, 1, 0]))

        formal_result = chargefw.calculate(molecule, method="formal", execution="full")
        self.assertEqual(formal_result.status, "success")

        assessment = chargefw.assess(
            molecule,
            method="qeq",
            parameter_set="QEq_original",
            execution="full",
        )
        self.assertFalse(assessment.plans)
        invalid_geometry = [
            issue
            for rejection in assessment.rejections
            for issue in rejection.issues
            if issue.kind == "invalid_geometry"
        ]
        self.assertEqual(len(invalid_geometry), 1)
        self.assertEqual(invalid_geometry[0].atom_index, 1)
        self.assertEqual(invalid_geometry[0].conformer_index, 0)
        self.assertEqual(
            invalid_geometry[0].message,
            "molecule 1: method 'qeq', conformer 1: atom 2 (H, formal charge 0) has non-finite "
            "coordinates",
        )

    def test_invalid_selection_requests_raise_value_error(self) -> None:
        with self.assertRaises(ValueError):
            chargefw.assess(water(), method="not-a-method")

        with self.assertRaisesRegex(ValueError, "parameter_set requires an explicit method"):
            chargefw.assess(water(), parameter_set="QEq_original")

        with self.assertRaises(ValueError):
            chargefw.assess(water(), method="eem", parameter_set="not-a-parameter-set")

        with self.assertRaises(ValueError):
            chargefw.assess(water(), method="eem", parameter_set="QEq_original")

    def test_results_outlive_calculation_inputs(self) -> None:
        def calculate_owned_result() -> tuple[chargefw.CalculationResult, chargefw.Plan]:
            molecule = chargefw.Molecule(
                [8, 1, 1],
                source_name="owned",
                record_id="owned-record",
                atom_ids=["O", "H1", "H2"],
            )
            assessment = assess_formal(molecule)
            plan = assessment.default_plan
            if plan is None:
                self.fail("assessment must produce a default plan")
            return chargefw.calculate(molecule, plan), plan

        result, plan = calculate_owned_result()
        gc.collect()
        self.assertEqual(plan.policy.mode, "full")
        self.assertEqual(result.status, "success")
        self.assertEqual(result.assignments[0].source.record_id, "owned-record")
        self.assertEqual(result.assignments[0].atom_ids, ("O", "H1", "H2"))
        np.testing.assert_array_equal(result.assignments[0].values, [0.0, 0.0, 0.0])

        values = result.assignments[0].values
        del result
        gc.collect()
        np.testing.assert_array_equal(values, [0.0, 0.0, 0.0])
        self.assertFalse(values.flags.writeable)

    def test_private_binding_rejects_invalid_policy_strings(self) -> None:
        molecule = water()
        with self.assertRaises(ValueError):
            _native_calculation._make_assessment(
                [molecule._native],
                [None],
                [("", 0, "")],
                [None],
                "",
                chargefw.calculation._default_parameter_catalog(),
                None,
                None,
                {},
                False,
                cast(Any, "invalid"),
                None,
                20_000,
                80_000,
                0,
            )

    def test_invalid_options_are_rejected_early(self) -> None:
        invalid_options = (
            (ValueError, {"execution": "fast"}),
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
        if result.plan is None:
            self.fail("successful calculation must report plan provenance")
        self.assertEqual(
            dict(result.plan.options),
            {"initial_charges": "zero", "iters": 2},
        )

        for options in (
            {"iters": 0},
            {"initial_charges": "unknown"},
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
        self.assertIsInstance(methods, Mapping)
        self.assertIn("eem", methods)
        self.assertEqual(methods.get("eem"), methods["eem"])
        self.assertIsNone(methods.get("not-a-method"))
        self.assertEqual(tuple(methods), tuple(method.id for method in methods.values()))
        self.assertEqual(methods["eem"].id, "eem")
        with self.assertRaisesRegex(KeyError, "unknown method ID"):
            methods["not-a-method"]
        with self.assertRaises(TypeError):
            cast(Any, methods)["unexpected"] = methods["eem"]
        with self.assertRaises(TypeError):
            del cast(Any, methods)["eem"]

        eem = methods["eem"]
        self.assertEqual(eem.notes, "")
        self.assertTrue(eem.requires_coordinates)
        self.assertEqual(eem.time_complexity, "O(n^3)")
        self.assertEqual(eem.memory_complexity, "O(n^2)")
        self.assertTrue(eem.supports_cutoff)
        self.assertTrue(eem.supports_cover)
        self.assertEqual(
            tuple(eem.parameter_sets),
            tuple(chargefw.parameter_sets.for_method("eem")),
        )
        peoe = methods["peoe"]
        self.assertIn("initial_charges=formal", peoe.notes)
        self.assertEqual(len(peoe.options), 2)
        self.assertIsInstance(peoe.options, Mapping)
        self.assertEqual(peoe.options["iters"].id, "iters")
        self.assertEqual(peoe.options["iters"].type, "integer")
        self.assertEqual(peoe.options["iters"].default, 6)
        self.assertEqual(peoe.options["iters"].minimum, 1)
        initial_charges = peoe.options["initial_charges"]
        self.assertEqual(initial_charges.type, "string")
        self.assertEqual(initial_charges.default, "zero")
        self.assertEqual(initial_charges.choices, ("zero", "formal"))
        qeq = methods["qeq"]
        overlap = qeq.options["overlap_term"]
        self.assertEqual(overlap.type, "string")
        self.assertIn("Ohno", overlap.choices)

        self.assertIsInstance(chargefw.parameter_sets, Mapping)
        parameter_set = next(iter(eem.parameter_sets.values()))
        self.assertEqual(chargefw.parameter_sets[parameter_set.id], parameter_set)
        result = chargefw.calculate(
            water(),
            method=eem,
            parameter_set=parameter_set,
            execution="full",
        )
        self.assertEqual(result.status, "success")
        executed = result.plan
        if executed is None:
            self.fail("successful calculation must report plan provenance")
        self.assertEqual(executed.parameter_set, parameter_set)


if __name__ == "__main__":
    unittest.main()
