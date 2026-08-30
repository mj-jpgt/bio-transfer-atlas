"""Competitive human-versus-microbial contig labelling."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import StrEnum
from itertools import groupby
from typing import Iterable

from hostbias.schemas import AlignmentRow, SchemaError


class ContigLabel(StrEnum):
    HUMAN = "human"
    AMBIGUOUS = "ambiguous"
    NON_HUMAN = "non_human"


@dataclass(frozen=True)
class LabelThresholds:
    min_aligned_bp: int = 500
    min_identity: float = 0.95
    min_query_coverage: float = 0.50
    min_mapq: float = 20
    min_human_score_ratio: float = 1.05

    def __post_init__(self) -> None:
        if self.min_aligned_bp < 1:
            raise ValueError("min_aligned_bp must be positive")
        if not 0 <= self.min_identity <= 1:
            raise ValueError("min_identity must be in [0, 1]")
        if not 0 <= self.min_query_coverage <= 1:
            raise ValueError("min_query_coverage must be in [0, 1]")
        if self.min_mapq < 0:
            raise ValueError("min_mapq must be non-negative")
        if self.min_human_score_ratio <= 1:
            raise ValueError("min_human_score_ratio must be greater than 1")


@dataclass(frozen=True)
class ContigCall:
    sample_id: str
    contig_id: str
    contig_length: int
    label: ContigLabel
    human_score: float | None
    gtdb_score: float | None
    human_to_gtdb_ratio: float | None

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        result["label"] = self.label.value
        return result


def _passes_human_thresholds(row: AlignmentRow, thresholds: LabelThresholds) -> bool:
    return (
        row.target_domain == "human"
        and row.aligned_bp >= thresholds.min_aligned_bp
        and row.identity >= thresholds.min_identity
        and row.query_coverage >= thresholds.min_query_coverage
        and row.mapq >= thresholds.min_mapq
    )


def label_contigs(
    rows: Iterable[AlignmentRow], thresholds: LabelThresholds | None = None
) -> list[ContigCall]:
    """Make one deterministic call per sample/contig from best domain hits.

    Human-like contigs whose competitive score is too close to GTDB are retained
    as ``ambiguous`` and excluded from all scientific denominators.
    """

    thresholds = thresholds or LabelThresholds()
    ordered = sorted(rows, key=lambda row: (row.sample_id, row.contig_id))
    calls: list[ContigCall] = []
    for (sample_id, contig_id), group in groupby(
        ordered, key=lambda row: (row.sample_id, row.contig_id)
    ):
        hits = list(group)
        lengths = {row.contig_length for row in hits}
        if len(lengths) != 1:
            raise SchemaError(
                f"inconsistent lengths for contig {(sample_id, contig_id)}"
            )
        human_hits = [row for row in hits if _passes_human_thresholds(row, thresholds)]
        gtdb_hits = [row for row in hits if row.target_domain == "gtdb"]
        best_human = max(human_hits, key=lambda row: row.alignment_score, default=None)
        best_gtdb = max(gtdb_hits, key=lambda row: row.alignment_score, default=None)
        if best_human is None:
            label = ContigLabel.NON_HUMAN
            ratio = None
        elif best_gtdb is None or best_gtdb.alignment_score == 0:
            label = ContigLabel.HUMAN
            ratio = None
        else:
            ratio = best_human.alignment_score / best_gtdb.alignment_score
            label = (
                ContigLabel.HUMAN
                if ratio >= thresholds.min_human_score_ratio
                else ContigLabel.AMBIGUOUS
            )
        calls.append(
            ContigCall(
                sample_id=sample_id,
                contig_id=contig_id,
                contig_length=lengths.pop(),
                label=label,
                human_score=None if best_human is None else best_human.alignment_score,
                gtdb_score=None if best_gtdb is None else best_gtdb.alignment_score,
                human_to_gtdb_ratio=ratio,
            )
        )
    return calls
