from collections.abc import Sequence

from .._payloads import AssessmentReportPayload, ExecutionPlanPayload, ExecutionResultPayload
from .._types import ChargeCorrection, Execution
from .core import _NativeMolecule
from .parameters import _NativeParameterCatalog

class _NativePlan:
    def report(self) -> ExecutionPlanPayload: ...
    def calculate(self, max_threads: int | None = None) -> ExecutionResultPayload: ...

class _NativeAssessment:
    def report(self) -> AssessmentReportPayload: ...
    def plans(self) -> list[_NativePlan]: ...
    def calculate_default(self) -> ExecutionResultPayload: ...

def _make_assessment(
    molecules: Sequence[_NativeMolecule],
    molecule_collection_name: str,
    catalog: _NativeParameterCatalog,
    method_id: str | None,
    parameter_set_id: str | None,
    method_options: dict[str, dict[str, bool | int | float | str]],
    permissive_types: bool,
    execution: Execution,
    radius: float | None,
    charge_correction: ChargeCorrection | None,
    cutoff_threshold: int | None,
    cover_threshold: int | None,
    max_threads: int,
) -> _NativeAssessment: ...
