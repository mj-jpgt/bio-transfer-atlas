from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from hostbias.data_manifest import ManifestError
from hostbias.reference_manifest import (
    build_balanced_panel,
    validate_reference_manifest,
    verify_downloads,
)


def catalog() -> list[dict[str, str]]:
    rows = []
    for group in ("AFR", "EUR"):
        for suffix in ("01", "02", "03"):
            donor = f"{group}{suffix}"
            rows.append(
                {
                    "donor_id": donor,
                    "population_group": group,
                    "assembly_url": f"https://example.test/{donor}.fa.gz",
                    "sha256": hashlib.sha256(donor.encode()).hexdigest(),
                }
            )
    return list(reversed(rows))


def build_panel() -> list[dict[str, str]]:
    return build_balanced_panel(
        catalog(),
        per_group=2,
        holdout_per_group=1,
        chm13_url="https://example.test/chm13.fa.gz",
        chm13_sha256=hashlib.sha256(b"chm13").hexdigest(),
    )


def test_build_balanced_panel_is_deterministic_and_disjoint() -> None:
    panel = build_panel()
    report = validate_reference_manifest(panel, per_group=2, holdout_per_group=1)
    assert report["population_groups"] == ["AFR", "EUR"]
    assert [row["donor_id"] for row in panel if row["population_group"] == "AFR"] == [
        "AFR01",
        "AFR02",
        "AFR03",
    ]
    assert [
        row["donor_id"] for row in panel if row["held_out_control"] == "true"
    ] == ["AFR03", "EUR03"]


def test_build_rejects_group_without_enough_holdout_donors() -> None:
    with pytest.raises(ManifestError, match="need 4 donors"):
        build_balanced_panel(
            catalog(),
            per_group=3,
            holdout_per_group=1,
            chm13_url="https://example.test/chm13.fa.gz",
            chm13_sha256=hashlib.sha256(b"chm13").hexdigest(),
        )


def test_verify_downloads_reports_ok_missing_and_mismatch(tmp_path: Path) -> None:
    panel = build_panel()[:3]
    first = tmp_path / panel[0]["local_filename"]
    first.write_bytes(b"chm13")
    second = tmp_path / panel[1]["local_filename"]
    second.write_bytes(b"wrong")
    report = verify_downloads(panel, tmp_path)
    assert report["valid"] is False
    assert [check["status"] for check in report["checks"]] == [
        "ok",
        "checksum_mismatch",
        "missing",
    ]

