"""Positive/negative control evaluation for competitive labelling."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

from hostbias.labeling import ContigCall, ContigLabel
from hostbias.schemas import ControlTruthRow, SchemaError, assert_unique


@dataclass(frozen=True)
class ControlResult:
    human_n: int
    microbial_n: int
    human_sensitivity: float
    human_sensitivity_bp: float
    microbial_false_positive_rate: float
    microbial_false_positive_bp_rate: float
    passed: bool

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def evaluate_controls(
    calls: Iterable[ContigCall],
    truth_rows: Iterable[ControlTruthRow],
    min_sensitivity: float = 0.95,
    max_false_positive_bp_rate: float = 0.001,
) -> ControlResult:
    """Require sensitivity by fragment and specificity by microbial base pairs."""

    calls = list(calls)
    truth_rows = list(truth_rows)
    assert_unique(calls, ("sample_id", "contig_id"))
    assert_unique(truth_rows, ("sample_id", "contig_id"))
    calls_by_key = {(row.sample_id, row.contig_id): row for row in calls}
    truth_keys = {(row.sample_id, row.contig_id) for row in truth_rows}
    if set(calls_by_key) != truth_keys:
        missing = truth_keys - set(calls_by_key)
        extra = set(calls_by_key) - truth_keys
        raise SchemaError(f"control call/truth key mismatch; missing={missing}, extra={extra}")
    humans = [row for row in truth_rows if row.truth == "human"]
    microbes = [row for row in truth_rows if row.truth == "microbial"]
    if not humans or not microbes:
        raise SchemaError("controls require both human and microbial truth")

    human_tp = [
        row
        for row in humans
        if calls_by_key[(row.sample_id, row.contig_id)].label == ContigLabel.HUMAN
    ]
    microbial_fp = [
        row
        for row in microbes
        if calls_by_key[(row.sample_id, row.contig_id)].label == ContigLabel.HUMAN
    ]
    sensitivity = len(human_tp) / len(humans)
    sensitivity_bp = sum(row.contig_length for row in human_tp) / sum(
        row.contig_length for row in humans
    )
    fp_rate = len(microbial_fp) / len(microbes)
    fp_bp_rate = sum(row.contig_length for row in microbial_fp) / sum(
        row.contig_length for row in microbes
    )
    return ControlResult(
        human_n=len(humans),
        microbial_n=len(microbes),
        human_sensitivity=sensitivity,
        human_sensitivity_bp=sensitivity_bp,
        microbial_false_positive_rate=fp_rate,
        microbial_false_positive_bp_rate=fp_bp_rate,
        passed=sensitivity >= min_sensitivity
        and fp_bp_rate <= max_false_positive_bp_rate,
    )
