"""Deterministic Gate A bootstrap, permutation, and influence analyses."""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np

from hostbias.endpoints import SampleEndpoint
from hostbias.schemas import SampleGroupRow, SchemaError, assert_unique


@dataclass(frozen=True)
class GroupedEndpoint:
    sample_id: str
    cohort: str
    p_count: float
    p_bp: float
    propagated_human_contig_count: int


@dataclass(frozen=True)
class MetricStatistics:
    metric: str
    tanzania_mean: float
    netherlands_mean: float
    difference: float
    ratio: float | None
    ratio_is_infinite: bool
    difference_ci_low: float
    difference_ci_high: float
    ratio_ci_low: float | None
    ratio_ci_high: float | None
    permutation_p_one_sided: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class AnalysisStatistics:
    p_count: MetricStatistics
    p_bp: MetricStatistics
    leave_one_out_preserves_p_count_direction: bool
    leave_one_out_min_difference: float
    max_tanzania_numerator_share: float
    tanzania_positive_samples: int
    group_sizes: dict[str, int]
    bootstrap_iterations: int
    permutation_iterations: int
    seed: int

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        return result


def attach_groups(
    endpoints: Iterable[SampleEndpoint], groups: Iterable[SampleGroupRow]
) -> list[GroupedEndpoint]:
    endpoints = list(endpoints)
    groups = list(groups)
    assert_unique(endpoints, ("sample_id",))
    assert_unique(groups, ("sample_id",))
    endpoint_ids = {row.sample_id for row in endpoints}
    group_by_sample = {row.sample_id: row.cohort for row in groups}
    if endpoint_ids != set(group_by_sample):
        raise SchemaError(
            "sample endpoint/group mismatch; "
            f"missing_groups={endpoint_ids - set(group_by_sample)}, "
            f"missing_endpoints={set(group_by_sample) - endpoint_ids}"
        )
    return [
        GroupedEndpoint(
            sample_id=row.sample_id,
            cohort=group_by_sample[row.sample_id],
            p_count=row.p_count,
            p_bp=row.p_bp,
            propagated_human_contig_count=row.propagated_human_contig_count,
        )
        for row in endpoints
    ]


def _finite_ratio(numerator: float, denominator: float) -> tuple[float | None, bool]:
    if denominator == 0:
        return (None, numerator > 0)
    return numerator / denominator, False


def _serializable_quantile(values: np.ndarray, q: float) -> float | None:
    ordered = np.sort(values)
    position = (len(ordered) - 1) * q
    lower_index = math.floor(position)
    upper_index = math.ceil(position)
    lower = float(ordered[lower_index])
    upper = float(ordered[upper_index])
    if not math.isfinite(lower) or not math.isfinite(upper):
        return None
    weight = position - lower_index
    return lower + weight * (upper - lower)


def _metric_statistics(
    tanzania: np.ndarray,
    netherlands: np.ndarray,
    metric: str,
    rng: np.random.Generator,
    bootstrap_iterations: int,
    permutation_iterations: int,
) -> MetricStatistics:
    t_mean = float(tanzania.mean())
    n_mean = float(netherlands.mean())
    difference = t_mean - n_mean
    ratio, ratio_is_infinite = _finite_ratio(t_mean, n_mean)

    t_indices = rng.integers(
        0, len(tanzania), size=(bootstrap_iterations, len(tanzania))
    )
    n_indices = rng.integers(
        0, len(netherlands), size=(bootstrap_iterations, len(netherlands))
    )
    boot_t = tanzania[t_indices].mean(axis=1)
    boot_n = netherlands[n_indices].mean(axis=1)
    boot_difference = boot_t - boot_n
    with np.errstate(divide="ignore", invalid="ignore"):
        boot_ratio = np.divide(
            boot_t,
            boot_n,
            out=np.full_like(boot_t, np.inf),
            where=boot_n != 0,
        )
        boot_ratio[(boot_t == 0) & (boot_n == 0)] = 1.0

    combined = np.concatenate((tanzania, netherlands))
    exceedances = 0
    for _ in range(permutation_iterations):
        permuted = rng.permutation(combined)
        permuted_difference = (
            permuted[: len(tanzania)].mean() - permuted[len(tanzania) :].mean()
        )
        if permuted_difference >= difference - 1e-15:
            exceedances += 1
    permutation_p = (exceedances + 1) / (permutation_iterations + 1)

    return MetricStatistics(
        metric=metric,
        tanzania_mean=t_mean,
        netherlands_mean=n_mean,
        difference=difference,
        ratio=ratio,
        ratio_is_infinite=ratio_is_infinite,
        difference_ci_low=float(np.quantile(boot_difference, 0.025)),
        difference_ci_high=float(np.quantile(boot_difference, 0.975)),
        ratio_ci_low=_serializable_quantile(boot_ratio, 0.025),
        ratio_ci_high=_serializable_quantile(boot_ratio, 0.975),
        permutation_p_one_sided=permutation_p,
    )


def analyze(
    rows: Iterable[GroupedEndpoint],
    bootstrap_iterations: int = 50_000,
    permutation_iterations: int = 100_000,
    seed: int = 20_260_729,
) -> AnalysisStatistics:
    """Calculate all primary and required secondary Gate A statistics."""

    rows = list(rows)
    if bootstrap_iterations < 100:
        raise ValueError("bootstrap_iterations must be at least 100")
    if permutation_iterations < 100:
        raise ValueError("permutation_iterations must be at least 100")
    assert_unique(rows, ("sample_id",))
    t_rows = [row for row in rows if row.cohort == "tanzania"]
    n_rows = [row for row in rows if row.cohort == "netherlands"]
    if not t_rows or not n_rows:
        raise SchemaError("both cohorts require at least one sample")
    unsupported = {row.cohort for row in rows} - {"tanzania", "netherlands"}
    if unsupported:
        raise SchemaError(f"unsupported cohorts: {unsupported}")
    for row in rows:
        if (
            not math.isfinite(row.p_count)
            or not math.isfinite(row.p_bp)
            or not 0 <= row.p_count <= 1
            or not 0 <= row.p_bp <= 1
        ):
            raise SchemaError(f"invalid propagation rate for sample {row.sample_id}")
        if row.propagated_human_contig_count < 0:
            raise SchemaError(
                f"negative propagated count for sample {row.sample_id}"
            )

    seed_sequence = np.random.SeedSequence(seed)
    rng_count, rng_bp = [
        np.random.default_rng(child) for child in seed_sequence.spawn(2)
    ]
    p_count = _metric_statistics(
        np.array([row.p_count for row in t_rows]),
        np.array([row.p_count for row in n_rows]),
        "p_count",
        rng_count,
        bootstrap_iterations,
        permutation_iterations,
    )
    p_bp = _metric_statistics(
        np.array([row.p_bp for row in t_rows]),
        np.array([row.p_bp for row in n_rows]),
        "p_bp",
        rng_bp,
        bootstrap_iterations,
        permutation_iterations,
    )

    loo_differences: list[float] = []
    for omitted in rows:
        kept_t = [row.p_count for row in t_rows if row.sample_id != omitted.sample_id]
        kept_n = [row.p_count for row in n_rows if row.sample_id != omitted.sample_id]
        if kept_t and kept_n:
            loo_differences.append(float(np.mean(kept_t) - np.mean(kept_n)))
    numerator_total = sum(row.propagated_human_contig_count for row in t_rows)
    max_share = (
        max(row.propagated_human_contig_count for row in t_rows) / numerator_total
        if numerator_total
        else 0.0
    )
    return AnalysisStatistics(
        p_count=p_count,
        p_bp=p_bp,
        leave_one_out_preserves_p_count_direction=all(
            value > 0 for value in loo_differences
        ),
        leave_one_out_min_difference=min(loo_differences),
        max_tanzania_numerator_share=max_share,
        tanzania_positive_samples=sum(
            row.propagated_human_contig_count > 0 for row in t_rows
        ),
        group_sizes={"tanzania": len(t_rows), "netherlands": len(n_rows)},
        bootstrap_iterations=bootstrap_iterations,
        permutation_iterations=permutation_iterations,
        seed=seed,
    )
