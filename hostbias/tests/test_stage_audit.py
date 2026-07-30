import gzip
import hashlib
import json
from pathlib import Path

import pytest

from hostbias.stage_audit import StageAuditError, audit_stage, inspect_fastq_pair


def _write_pair(
    root: Path,
    stem: str,
    names: list[str],
    *,
    length: int = 4,
) -> tuple[Path, Path]:
    paths = root / f"{stem}_R1.fastq.gz", root / f"{stem}_R2.fastq.gz"
    for mate, path in enumerate(paths, start=1):
        with gzip.open(path, "wt", encoding="ascii") as handle:
            for name in names:
                sequence = ("ACGT" * ((length + 3) // 4))[:length]
                handle.write(f"@{name}/{mate}\n{sequence}\n+\n{'I' * length}\n")
    return paths


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_stage_audit_is_aggregate_safe_and_exact(tmp_path: Path) -> None:
    normalized = _write_pair(tmp_path, "normalized", ["read-a", "read-b", "read-c"])
    source = _write_pair(tmp_path, "source", ["read-a", "read-b"])
    strict = _write_pair(tmp_path, "strict", ["read-a"])

    report = audit_stage(
        sample_id="SRR1",
        normalized_r1=normalized[0],
        normalized_r2=normalized[1],
        source_r1=source[0],
        source_r2=source[1],
        strict_r1=strict[0],
        strict_r2=strict[1],
        expected_r1_sha256=_sha256(normalized[0]),
        expected_r2_sha256=_sha256(normalized[1]),
        expected_r1_bytes=normalized[0].stat().st_size,
        expected_r2_bytes=normalized[1].stat().st_size,
        expected_pairs=3,
        expected_length=4,
    )

    assert report["status"] == "PASS"
    assert report["normalized"]["pairs"] == 3
    assert report["normalized"]["matches_prefilter_baseline"] is True
    assert report["grch38_filter"]["source"]["removed_pairs"] == 1
    assert report["grch38_filter"]["strict"]["removed_pairs"] == 2
    assert report["grch38_filter"]["source"]["removed_fraction"] == pytest.approx(1 / 3)
    assert report["privacy"]["read_identifiers_recorded"] is False
    serialized = json.dumps(report)
    assert str(tmp_path) not in serialized
    assert "read-a" not in serialized
    assert "ACGT" not in serialized


def test_stage_audit_rejects_changed_normalized_input(tmp_path: Path) -> None:
    normalized = _write_pair(tmp_path, "normalized", ["read-a"])
    source = _write_pair(tmp_path, "source", ["read-a"])
    strict = _write_pair(tmp_path, "strict", ["read-a"])

    with pytest.raises(StageAuditError, match="pre-filter baseline"):
        audit_stage(
            sample_id="SRR1",
            normalized_r1=normalized[0],
            normalized_r2=normalized[1],
            source_r1=source[0],
            source_r2=source[1],
            strict_r1=strict[0],
            strict_r2=strict[1],
            expected_r1_sha256="0" * 64,
            expected_r2_sha256=_sha256(normalized[1]),
            expected_r1_bytes=normalized[0].stat().st_size,
            expected_r2_bytes=normalized[1].stat().st_size,
            expected_pairs=1,
            expected_length=4,
        )


def test_pair_audit_rejects_desynchronization_without_leaking_names(
    tmp_path: Path,
) -> None:
    r1, _ = _write_pair(tmp_path, "first", ["sensitive-a"])
    _, r2 = _write_pair(tmp_path, "second", ["sensitive-b"])

    with pytest.raises(StageAuditError) as error:
        inspect_fastq_pair(r1, r2, expected_pairs=1, expected_length=4)

    assert str(error.value) == "mate identifiers differ at pair 1"
    assert "sensitive" not in str(error.value)
