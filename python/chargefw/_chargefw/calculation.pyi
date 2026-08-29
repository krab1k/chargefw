from collections.abc import Sequence

from .._payloads import AssessmentReportPayload, ExecutionResultPayload
from .._types import ChargeCorrection, Execution
from .core import _NativeMolecule
from .parameters import _NativeParameterCatalog

class _NativeAssessment:
    def report(self) -> AssessmentReportPayload: ...
    def calculate(self) -> ExecutionResultPayload: ...

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
