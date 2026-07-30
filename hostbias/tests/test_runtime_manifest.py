from __future__ import annotations

import csv
from pathlib import Path

import pytest

from hostbias.data_manifest import ManifestError
from hostbias.runtime_manifest import build_runtime_rows


PROJECT = Path(__file__).parents[1]


def frozen_rows() -> list[dict[str, str]]:
    with (PROJECT / "config" / "stage0_samples.tsv").open(
        encoding="utf-8", newline=""
    ) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def snapshots(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = {}
    for frozen in rows:
        run = frozen["run_accession"]
        result.setdefault(frozen["arm"], []).append(
            {
                "run_accession": run,
                "study_accession": frozen["bioproject"],
                "sample_accession": frozen["sample_accession"],
                "library_layout": frozen["library_layout"],
                "instrument_platform": frozen["instrument_platform"],
                "library_strategy": frozen["library_strategy"],
                "base_count": frozen["base_count"],
                "read_count": frozen["read_count"],
                "fastq_ftp": f"ftp.example/{run}_1.fastq.gz;ftp.example/{run}_2.fastq.gz",
                "fastq_md5": "1" * 32 + ";" + "2" * 32,
                "fastq_bytes": "100;200",
            }
        )
    for arm in result:
        result[arm].sort(key=lambda row: row["run_accession"])
    return result


def test_sentinel_scope_preserves_frozen_arm_and_rank_order() -> None:
    frozen = frozen_rows()
    rows = build_runtime_rows(frozen, snapshots(frozen), scope="sentinel")

    assert len(rows) == 6
    assert [row["accession"] for row in rows] == [
        row["run_accession"]
        for row in frozen
        if row["role"] == "primary" and int(row["rank"]) <= 3
    ]
    assert {row["status"] for row in rows} == {"primary"}
    assert rows[0]["fastq1_url"].startswith("ftp://")


def test_primary_scope_has_exactly_twenty_runs_per_arm() -> None:
    frozen = frozen_rows()
    rows = build_runtime_rows(frozen, snapshots(frozen), scope="primary")

    assert len(rows) == 40
    assert sum(row["cohort"] == "tanzania" for row in rows) == 20
    assert sum(row["cohort"] == "netherlands" for row in rows) == 20


def test_duplicate_snapshot_accession_is_rejected() -> None:
    frozen = frozen_rows()
    metadata = snapshots(frozen)
    metadata["tanzania"].append(dict(metadata["tanzania"][0]))

    with pytest.raises(ManifestError, match="duplicate run accessions"):
        build_runtime_rows(frozen, metadata, scope="sentinel")


def test_snapshot_metadata_drift_is_rejected() -> None:
    frozen = frozen_rows()
    metadata = snapshots(frozen)
    selected = next(
        row
        for row in metadata["tanzania"]
        if row["run_accession"] == frozen[0]["run_accession"]
    )
    selected["sample_accession"] = "SAMN_CHANGED"

    with pytest.raises(ManifestError, match="differs from ENA snapshot"):
        build_runtime_rows(frozen, metadata, scope="sentinel")


def test_noncanonical_snapshot_order_is_rejected() -> None:
    frozen = frozen_rows()
    metadata = snapshots(frozen)
    metadata["netherlands"][0], metadata["netherlands"][1] = (
        metadata["netherlands"][1],
        metadata["netherlands"][0],
    )

    with pytest.raises(ManifestError, match="canonical accession order"):
        build_runtime_rows(frozen, metadata, scope="sentinel")


def test_missing_size_or_wrong_file_count_is_rejected() -> None:
    frozen = frozen_rows()
    metadata = snapshots(frozen)
    selected = next(
        row
        for row in metadata["tanzania"]
        if row["run_accession"] == frozen[0]["run_accession"]
    )
    selected["fastq_bytes"] = "100"

    with pytest.raises(ManifestError, match="exactly two values"):
        build_runtime_rows(frozen, metadata, scope="sentinel")
