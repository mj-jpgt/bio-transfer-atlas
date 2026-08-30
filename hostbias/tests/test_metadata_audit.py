from __future__ import annotations

import hashlib
from pathlib import Path

from hostbias.data_manifest import MANIFEST_FIELDS, canonical_tsv, select_arm
from hostbias.metadata_audit import audit_live_metadata


def row(accession: str, sample: str, bases: int) -> dict[str, str]:
    return {
        "run_accession": accession,
        "study_accession": "ERPTEST",
        "sample_accession": sample,
        "library_layout": "PAIRED",
        "instrument_platform": "ILLUMINA",
        "library_strategy": "WGS",
        "base_count": str(bases),
        "read_count": str(bases // 100),
        "fastq_ftp": f"ftp.test/{accession}.fastq.gz",
        "fastq_md5": hashlib.md5(accession.encode()).hexdigest(),
    }


def catalogs() -> dict[str, list[dict[str, str]]]:
    return {
        "PRJA": [row("A1", "AS1", 500), row("A2", "AS2", 400), row("A3", "AS3", 300)],
        "PRJB": [row("B1", "BS1", 700), row("B2", "BS2", 600), row("B3", "BS3", 200)],
    }


def write_frozen(path: Path, source: dict[str, list[dict[str, str]]]) -> None:
    selected = []
    for arm, project in (("a", "PRJA"), ("b", "PRJB")):
        selected.extend(
            select_arm(
                source[project],
                arm=arm,
                project=project,
                primary=2,
                reserves=1,
            )
        )
    path.write_bytes(canonical_tsv(selected, MANIFEST_FIELDS))


def test_live_audit_matches_frozen_manifest(tmp_path: Path) -> None:
    source = catalogs()
    frozen = tmp_path / "frozen.tsv"
    write_frozen(frozen, source)
    report = audit_live_metadata(
        frozen,
        [("a", "PRJA"), ("b", "PRJB")],
        primary=2,
        reserves=1,
        fetcher=lambda project: source[project],
        checked_at_utc="2026-07-30T00:00:00+00:00",
    )
    assert report["valid"] is True
    assert report["manifest_exact_match"] is True
    assert report["frozen_manifest_sha256"] == report["live_manifest_sha256"]
    assert report["arms"]["a"]["response_rows"] == 3


def test_live_audit_detects_selection_drift(tmp_path: Path) -> None:
    source = catalogs()
    frozen = tmp_path / "frozen.tsv"
    write_frozen(frozen, source)
    source["PRJA"][2]["base_count"] = "900"
    report = audit_live_metadata(
        frozen,
        [("a", "PRJA"), ("b", "PRJB")],
        primary=2,
        reserves=1,
        fetcher=lambda project: source[project],
        checked_at_utc="2026-07-30T00:00:00+00:00",
    )
    assert report["valid"] is False
    assert report["manifest_exact_match"] is False
    assert report["arms"]["a"]["selection_matches_frozen"] is False
