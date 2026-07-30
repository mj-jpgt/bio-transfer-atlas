"""Machine-readable, deterministic Gate A decision."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from hostbias.controls import ControlResult
from hostbias.schemas import SensitivityRow
from hostbias.statistics import AnalysisStatistics


@dataclass(frozen=True)
class GateThresholds:
    samples_per_cohort: int = 20
    min_tanzania_p_count: float = 0.01
    min_ratio: float = 1.5
    max_permutation_p: float = 0.05
    min_positive_tanzania_samples: int = 5
    max_sample_numerator_share: float = 0.25
    required_sensitivity_analyses: tuple[str, ...] = (
        "strict_pair",
        "identity_0.90",
        "identity_0.95",
        "identity_0.98",
    )


@dataclass(frozen=True)
class Criterion:
    name: str
    passed: bool
    observed: object
    required: str


@dataclass(frozen=True)
class GateVerdict:
    schema_version: str
    status: str
    pass_gate_a: bool | None
    first_failed_criterion: str | None
    criteria: tuple[Criterion, ...]
    statistics: AnalysisStatistics
    controls: ControlResult

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def make_verdict(
    statistics: AnalysisStatistics,
    controls: ControlResult,
    sensitivities: Iterable[SensitivityRow],
    thresholds: GateThresholds | None = None,
) -> GateVerdict:
    """Return PASS/FAIL only after operational prerequisites are satisfied."""

    thresholds = thresholds or GateThresholds()
    sensitivities = list(sensitivities)
    sensitivity_by_id = {
        (row.analysis_id, row.metric): row for row in sensitivities
    }
    required_keys = {
        (analysis_id, "p_count")
        for analysis_id in thresholds.required_sensitivity_analyses
    }
    missing_sensitivities = sorted(required_keys - set(sensitivity_by_id))
    operational_criteria = [
        Criterion(
            "controls_pass",
            controls.passed,
            controls.to_dict(),
            "sensitivity >= 0.95 and microbial false-positive bp rate <= 0.001",
        ),
        Criterion(
            "complete_groups",
            statistics.group_sizes
            == {
                "tanzania": thresholds.samples_per_cohort,
                "netherlands": thresholds.samples_per_cohort,
            },
            statistics.group_sizes,
            f"exactly {thresholds.samples_per_cohort} valid samples per cohort",
        ),
        Criterion(
            "complete_sensitivity_matrix",
            not missing_sensitivities,
            missing_sensitivities,
            f"p_count results for {list(thresholds.required_sensitivity_analyses)}",
        ),
    ]
    failed_operational = next(
        (criterion for criterion in operational_criteria if not criterion.passed), None
    )
    if failed_operational is not None:
        return GateVerdict(
            schema_version="1.0",
            status="OPERATIONAL_FAILURE",
            pass_gate_a=None,
            first_failed_criterion=failed_operational.name,
            criteria=tuple(operational_criteria),
            statistics=statistics,
            controls=controls,
        )

    primary = statistics.p_count
    ratio_value = float("inf") if primary.ratio_is_infinite else primary.ratio
    ratio_ci_low = (
        float("inf")
        if primary.ratio_ci_low is None and primary.ratio_is_infinite
        else primary.ratio_ci_low
    )
    sensitivity_direction = all(
        sensitivity_by_id[(analysis_id, "p_count")].tanzania_mean
        > sensitivity_by_id[(analysis_id, "p_count")].netherlands_mean
        for analysis_id in thresholds.required_sensitivity_analyses
    )
    criteria = operational_criteria + [
        Criterion(
            "minimum_tanzania_propagation",
            primary.tanzania_mean >= thresholds.min_tanzania_p_count,
            primary.tanzania_mean,
            f">= {thresholds.min_tanzania_p_count}",
        ),
        Criterion(
            "minimum_mean_ratio",
            ratio_value is not None and ratio_value >= thresholds.min_ratio,
            primary.ratio if not primary.ratio_is_infinite else "infinity",
            f">= {thresholds.min_ratio}",
        ),
        Criterion(
            "ratio_ci_excludes_one",
            ratio_ci_low is not None and ratio_ci_low > 1,
            primary.ratio_ci_low if primary.ratio_ci_low is not None else "infinity",
            "> 1",
        ),
        Criterion(
            "difference_ci_excludes_zero",
            primary.difference_ci_low > 0,
            primary.difference_ci_low,
            "> 0",
        ),
        Criterion(
            "permutation_significant",
            primary.permutation_p_one_sided < thresholds.max_permutation_p,
            primary.permutation_p_one_sided,
            f"< {thresholds.max_permutation_p}",
        ),
        Criterion(
            "minimum_positive_tanzania_samples",
            statistics.tanzania_positive_samples
            >= thresholds.min_positive_tanzania_samples,
            statistics.tanzania_positive_samples,
            f">= {thresholds.min_positive_tanzania_samples}",
        ),
        Criterion(
            "p_bp_preserves_direction",
            statistics.p_bp.difference > 0,
            statistics.p_bp.difference,
            "> 0",
        ),
        Criterion(
            "sensitivity_matrix_preserves_direction",
            sensitivity_direction,
            {
                key[0]: {
                    "tanzania": row.tanzania_mean,
                    "netherlands": row.netherlands_mean,
                }
                for key, row in sensitivity_by_id.items()
                if key[1] == "p_count"
            },
            "Tanzania mean > Netherlands mean in every required analysis",
        ),
        Criterion(
            "leave_one_out_preserves_direction",
            statistics.leave_one_out_preserves_p_count_direction,
            statistics.leave_one_out_min_difference,
            "all leave-one-sample-out differences > 0",
        ),
        Criterion(
            "no_dominant_sample",
            statistics.max_tanzania_numerator_share
            <= thresholds.max_sample_numerator_share,
            statistics.max_tanzania_numerator_share,
            f"<= {thresholds.max_sample_numerator_share}",
        ),
    ]
    first_failure = next(
        (criterion.name for criterion in criteria if not criterion.passed), None
    )
    return GateVerdict(
        schema_version="1.0",
        status="PASS" if first_failure is None else "FAIL",
        pass_gate_a=first_failure is None,
        first_failed_criterion=first_failure,
        criteria=tuple(criteria),
        statistics=statistics,
        controls=controls,
    )


def write_verdict(verdict: GateVerdict, output_dir: str | Path) -> tuple[Path, Path]:
    """Write strict JSON and a concise human-readable sibling report."""

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "GATE_A_VERDICT.json"
    markdown_path = output_dir / "GATE_A_VERDICT.md"
    json_path.write_text(
        json.dumps(verdict.to_dict(), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    lines = [
        "# Gate A verdict",
        "",
        f"**Status:** {verdict.status}",
        "",
        f"First failed criterion: `{verdict.first_failed_criterion or 'none'}`",
        "",
        "| Criterion | Pass | Observed | Required |",
        "|---|---:|---|---|",
    ]
    lines.extend(
        f"| {row.name} | {'yes' if row.passed else 'no'} | "
        f"`{json.dumps(row.observed, sort_keys=True)}` | {row.required} |"
        for row in verdict.criteria
    )
    markdown_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path, markdown_path
