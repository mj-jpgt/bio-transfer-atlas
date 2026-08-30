"""Aggregate-only closeout audit for atomically published FASTQ pairs."""

from __future__ import annotations

import csv
import hashlib
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path


REQUIRED_FIELDS = {
    "sample_id",
    "fastq1_md5",
    "fastq1_bytes",
    "fastq2_md5",
    "fastq2_bytes",
}


class FetchAuditError(ValueError):
    """Raised when the published acquisition does not match its manifest."""


def _md5_stable(path: Path) -> tuple[int, str, bool]:
    before = path.stat()
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    after = path.stat()
    stable = (
        before.st_size,
        before.st_mtime_ns,
        before.st_ino,
    ) == (
        after.st_size,
        after.st_mtime_ns,
        after.st_ino,
    )
    return before.st_size, digest.hexdigest(), stable


def audit_fetch(
    *,
    manifest: Path,
    raw_root: Path,
    expected_pairs: int = 40,
    threads: int = 4,
) -> dict[str, object]:
    """Verify every expected final mate while excluding partial files."""

    manifest_bytes = manifest.read_bytes()
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or REQUIRED_FIELDS.difference(reader.fieldnames):
            raise FetchAuditError("fetch manifest is missing required fields")
        rows = list(reader)

    if len(rows) != expected_pairs:
        raise FetchAuditError(
            f"expected {expected_pairs} manifest pairs, observed {len(rows)}"
        )
    sample_ids = [row["sample_id"] for row in rows]
    if len(set(sample_ids)) != len(sample_ids):
        raise FetchAuditError("fetch manifest contains duplicate sample identifiers")

    expected: dict[Path, tuple[int, str]] = {}
    for row in rows:
        for mate in (1, 2):
            try:
                byte_count = int(row[f"fastq{mate}_bytes"])
            except ValueError as error:
                raise FetchAuditError("manifest contains a non-integer byte count") from error
            md5 = row[f"fastq{mate}_md5"].lower()
            if byte_count <= 0 or len(md5) != 32:
                raise FetchAuditError("manifest contains an invalid size or MD5")
            expected[raw_root / f"{row['sample_id']}_R{mate}.fastq.gz"] = (
                byte_count,
                md5,
            )

    actual_final = set(raw_root.glob("*.fastq.gz"))
    missing = set(expected).difference(actual_final)
    unexpected = actual_final.difference(expected)
    if missing or unexpected:
        raise FetchAuditError(
            "atomic final-file set mismatch: "
            f"missing={len(missing)}, unexpected={len(unexpected)}"
        )

    with ThreadPoolExecutor(max_workers=threads) as executor:
        observed = dict(zip(expected, executor.map(_md5_stable, expected), strict=True))

    size_matches = 0
    md5_matches = 0
    stable_files = 0
    total_bytes = 0
    for path, (expected_bytes, expected_md5) in expected.items():
        observed_bytes, observed_md5, stable = observed[path]
        total_bytes += observed_bytes
        size_matches += observed_bytes == expected_bytes
        md5_matches += observed_md5 == expected_md5
        stable_files += stable

    expected_files = expected_pairs * 2
    if (
        size_matches != expected_files
        or md5_matches != expected_files
        or stable_files != expected_files
    ):
        raise FetchAuditError(
            "published-file verification failed: "
            f"size_matches={size_matches}, md5_matches={md5_matches}, "
            f"stable_files={stable_files}, expected={expected_files}"
        )

    partial_files = sum(
        path.is_file() and ".partial" in path.name for path in raw_root.iterdir()
    )
    expected_bytes = sum(value[0] for value in expected.values())
    if total_bytes != expected_bytes:
        raise FetchAuditError("aggregate published byte count differs from manifest")

    return {
        "schema_version": 1,
        "checkpoint": "primary_fastq_acquisition_closeout",
        "captured_at_utc": datetime.now(UTC).isoformat(),
        "status": "PASS",
        "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "expected": {
            "pairs": expected_pairs,
            "mate_files": expected_files,
            "bytes": expected_bytes,
        },
        "observed": {
            "complete_pairs": expected_pairs,
            "final_mate_files": expected_files,
            "bytes": total_bytes,
            "size_matches": size_matches,
            "ena_md5_matches": md5_matches,
            "stable_during_audit": stable_files,
            "unexpected_final_files": 0,
            "partial_files_excluded": partial_files,
        },
        "privacy": {
            "aggregate_evidence_only": True,
            "absolute_paths_recorded": False,
            "sample_accessions_recorded": False,
            "filenames_recorded": False,
            "sequence_or_quality_recorded": False,
        },
    }
