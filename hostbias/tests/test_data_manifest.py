from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import pytest

from hostbias.data_manifest import (
    ENA_FIELDS,
    MANIFEST_FIELDS,
    ManifestError,
    select_arm,
    validate_manifest,
    write_snapshot,
)


def ena_row(
    accession: str,
    sample: str,
    bases: int,
    *,
    layout: str = "PAIRED",
    platform: str = "ILLUMINA",
    strategy: str = "WGS",
) -> dict[str, str]:
    return {
        "run_accession": accession,
        "study_accession": "PRJTEST",
        "sample_accession": sample,
        "library_layout": layout,
        "instrument_platform": platform,
        "library_strategy": strategy,
        "base_count": str(bases),
        "read_count": str(bases // 100),
        "fastq_ftp": f"ftp.test/{accession}_1.fastq.gz;ftp.test/{accession}_2.fastq.gz",
        "fastq_md5": "a;b",
    }


def test_select_arm_filters_deduplicates_and_breaks_ties_by_accession() -> None:
    rows = [
        ena_row("RUN_B", "S2", 500),
        ena_row("RUN_A", "S1", 500),
        ena_row("RUN_DUP", "S1", 400),
        ena_row("RUN_C", "S3", 300),
        ena_row("RUN_SINGLE", "S4", 900, layout="SINGLE"),
    ]
    selected = select_arm(
        rows, arm="test", project="PRJTEST", primary=2, reserves=1
    )
    assert [row["run_accession"] for row in selected] == ["RUN_A", "RUN_B", "RUN_C"]
    assert [row["role"] for row in selected] == ["primary", "primary", "reserve"]
    assert [row["rank"] for row in selected] == [1, 2, 3]


def test_select_arm_rejects_too_few_eligible_samples() -> None:
    with pytest.raises(ManifestError, match="need 3 eligible"):
        select_arm(
            [ena_row("RUN_A", "S1", 500)],
            arm="test",
            project="PRJTEST",
            primary=2,
            reserves=1,
        )


def test_snapshot_is_canonical_and_sidecar_matches(tmp_path: Path) -> None:
    output = tmp_path / "snapshot.tsv"
    rows = [ena_row("RUN_A", "S1", 500)]
    digest = write_snapshot(output, rows, ENA_FIELDS)
    assert digest == hashlib.sha256(output.read_bytes()).hexdigest()
    assert output.read_bytes().endswith(b"\n")
    assert output.with_suffix(".tsv.sha256").read_text(encoding="ascii") == (
        f"{digest}  snapshot.tsv\n"
    )


def test_validate_manifest_rejects_outcome_driven_reordering() -> None:
    selected = select_arm(
        [
            ena_row("RUN_A", "S1", 500),
            ena_row("RUN_B", "S2", 400),
            ena_row("RUN_C", "S3", 300),
        ],
        arm="test",
        project="PRJTEST",
        primary=2,
        reserves=1,
    )
    rows = [{key: str(row[key]) for key in MANIFEST_FIELDS} for row in selected]
    rows[0], rows[1] = rows[1], rows[0]
    with pytest.raises(ManifestError, match="ranks must be ordered"):
        validate_manifest(rows, primary=2, reserves=1)


def test_frozen_stage0_manifest_has_twenty_plus_ten_per_arm() -> None:
    path = Path(__file__).parents[1] / "config" / "stage0_samples.tsv"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    report = validate_manifest(rows)
    assert report["total_rows"] == 60
    assert set(report["arms"]) == {"tanzania", "netherlands"}

