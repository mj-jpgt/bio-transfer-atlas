import hashlib
import json
from pathlib import Path

import pytest

from hostbias.fetch_audit import FetchAuditError, audit_fetch


HEADER = (
    "sample_id\tfastq1_md5\tfastq1_bytes\tfastq2_md5\tfastq2_bytes\n"
)


def _prepare_pair(root: Path, sample: str, first: bytes, second: bytes) -> str:
    values: list[str] = [sample]
    for mate, payload in enumerate((first, second), start=1):
        path = root / f"{sample}_R{mate}.fastq.gz"
        path.write_bytes(payload)
        values.extend(
            [
                hashlib.md5(payload, usedforsecurity=False).hexdigest(),
                str(len(payload)),
            ]
        )
    return "\t".join(values)


def test_fetch_audit_verifies_all_finals_and_excludes_partials(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    rows = [
        _prepare_pair(raw, "SRR1", b"first-1", b"first-2"),
        _prepare_pair(raw, "SRR2", b"second-1", b"second-2"),
    ]
    (raw / ".SRR3_R1.fastq.gz.partial").write_bytes(b"incomplete")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")

    report = audit_fetch(
        manifest=manifest,
        raw_root=raw,
        expected_pairs=2,
        threads=2,
    )

    assert report["status"] == "PASS"
    assert report["observed"]["complete_pairs"] == 2
    assert report["observed"]["final_mate_files"] == 4
    assert report["observed"]["ena_md5_matches"] == 4
    assert report["observed"]["partial_files_excluded"] == 1
    serialized = json.dumps(report)
    assert str(tmp_path) not in serialized
    assert "SRR1" not in serialized
    assert "first-1" not in serialized


def test_fetch_audit_rejects_incomplete_pair_set(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    row = _prepare_pair(raw, "SRR1", b"first", b"second")
    (raw / "SRR1_R2.fastq.gz").unlink()
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(HEADER + row + "\n", encoding="utf-8")

    with pytest.raises(FetchAuditError, match="missing=1, unexpected=0"):
        audit_fetch(manifest=manifest, raw_root=raw, expected_pairs=1)


def test_fetch_audit_rejects_md5_mismatch_without_naming_file(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    row = _prepare_pair(raw, "SRR-sensitive", b"first", b"second")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        HEADER
        + row.replace(
            hashlib.md5(b"first", usedforsecurity=False).hexdigest(),
            "0" * 32,
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(FetchAuditError) as error:
        audit_fetch(manifest=manifest, raw_root=raw, expected_pairs=1)

    assert "md5_matches=1" in str(error.value)
    assert "sensitive" not in str(error.value)
