from __future__ import annotations

import pytest

from hostbias.data_manifest import ManifestError
from hostbias.sentinel import evaluate_sentinels


def manifest_rows() -> list[dict[str, str]]:
    rows = []
    for arm in ("a", "b"):
        for rank in range(1, 4):
            rows.append(
                {
                    "run_accession": f"{arm}{rank}",
                    "arm": arm,
                    "role": "primary",
                    "rank": str(rank),
                }
            )
    return rows


def metric(run: str, arm: str, mapped: int, **overrides: str) -> dict[str, str]:
    row = {
        "run_accession": run,
        "arm": arm,
        "streamed_spots": "1000000",
        "grch38_mapped_pairs": str(mapped),
        "metadata_ok": "true",
        "checksum_ok": "true",
    }
    row.update(overrides)
    return row


def passing_metrics() -> list[dict[str, str]]:
    return [
        metric("a1", "a", 100),
        metric("a2", "a", 101),
        metric("a3", "a", 99),
        metric("b1", "b", 100),
        metric("b2", "b", 100),
        metric("b3", "b", 0),
    ]


def test_exact_boundary_and_two_of_three_pass() -> None:
    report = evaluate_sentinels(manifest_rows(), passing_metrics())
    assert report["eligible"] is True
    assert report["arms"]["a"]["passing_host_runs"] == 2
    assert report["arms"]["a"]["runs"][0]["mapped_pairs_per_million"] == 100


def test_any_provenance_failure_makes_arm_ineligible() -> None:
    metrics = passing_metrics()
    metrics[2]["checksum_ok"] = "false"
    report = evaluate_sentinels(manifest_rows(), metrics)
    assert report["eligible"] is False
    assert report["arms"]["a"]["all_metadata_and_checksums_ok"] is False


def test_missing_frozen_sentinel_is_operational_error() -> None:
    with pytest.raises(ManifestError, match="missing sentinel metrics"):
        evaluate_sentinels(manifest_rows(), passing_metrics()[:-1])

