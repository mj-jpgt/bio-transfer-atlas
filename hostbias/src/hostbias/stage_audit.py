"""Privacy-safe audits for normalized and host-filtered paired FASTQ stages."""

from __future__ import annotations

import gzip
import hashlib
import re
from contextlib import ExitStack
from datetime import UTC, datetime
from pathlib import Path
from typing import TextIO


SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


class StageAuditError(ValueError):
    """Raised when a production-stage file contract is violated."""


def _open_fastq(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii")
    return path.open("r", encoding="ascii")


def _canonical_name(header: str, pair_number: int) -> str:
    if not header.startswith("@"):
        raise StageAuditError(f"invalid FASTQ header at pair {pair_number}")
    name = header[1:].split()[0]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    if not name:
        raise StageAuditError(f"empty FASTQ identifier at pair {pair_number}")
    return name


def _read_record(handle: TextIO, pair_number: int) -> tuple[str, int] | None:
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline().rstrip("\r\n")
    separator = handle.readline()
    quality = handle.readline().rstrip("\r\n")
    if not separator.startswith("+") or len(sequence) != len(quality):
        raise StageAuditError(f"invalid FASTQ record at pair {pair_number}")
    return _canonical_name(header, pair_number), len(sequence)


def inspect_fastq_pair(
    r1: Path,
    r2: Path,
    *,
    expected_pairs: int | None = None,
    expected_length: int | None = None,
) -> dict[str, int | bool]:
    """Stream a pair without retaining identifiers or sequence."""

    count = 0
    observed_length: int | None = None
    with ExitStack() as stack:
        first = stack.enter_context(_open_fastq(r1))
        second = stack.enter_context(_open_fastq(r2))
        while True:
            pair_number = count + 1
            record1 = _read_record(first, pair_number)
            record2 = _read_record(second, pair_number)
            if record1 is None and record2 is None:
                break
            if record1 is None or record2 is None:
                raise StageAuditError(
                    f"mate files contain different record counts at pair {pair_number}"
                )
            if record1[0] != record2[0]:
                raise StageAuditError(f"mate identifiers differ at pair {pair_number}")
            if record1[1] != record2[1]:
                raise StageAuditError(f"mate lengths differ at pair {pair_number}")
            if expected_length is not None and record1[1] != expected_length:
                raise StageAuditError(
                    f"unexpected read length at pair {pair_number}: "
                    f"expected {expected_length}, observed {record1[1]}"
                )
            if observed_length is None:
                observed_length = record1[1]
            elif observed_length != record1[1]:
                raise StageAuditError(f"non-uniform read length at pair {pair_number}")
            count += 1

    if expected_pairs is not None and count != expected_pairs:
        raise StageAuditError(f"expected {expected_pairs} pairs, observed {count}")
    return {
        "pairs": count,
        "read_length": observed_length if observed_length is not None else 0,
        "synchronized": True,
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _stable_file_state(path: Path) -> tuple[int, int, int]:
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns, stat.st_ino


def _removal_metrics(normalized_pairs: int, retained_pairs: int) -> dict[str, int | float]:
    if retained_pairs > normalized_pairs:
        raise StageAuditError(
            f"filtered pair count {retained_pairs} exceeds normalized count {normalized_pairs}"
        )
    removed_pairs = normalized_pairs - retained_pairs
    return {
        "retained_pairs": retained_pairs,
        "removed_pairs": removed_pairs,
        "retained_fraction": retained_pairs / normalized_pairs,
        "removed_fraction": removed_pairs / normalized_pairs,
    }


def audit_stage(
    *,
    sample_id: str,
    normalized_r1: Path,
    normalized_r2: Path,
    source_r1: Path,
    source_r2: Path,
    strict_r1: Path,
    strict_r2: Path,
    expected_r1_sha256: str,
    expected_r2_sha256: str,
    expected_r1_bytes: int,
    expected_r2_bytes: int,
    expected_pairs: int = 8_000_000,
    expected_length: int = 100,
) -> dict[str, object]:
    """Audit one completed normalization/GRCh38-filter stage."""

    expected_hashes = (expected_r1_sha256.lower(), expected_r2_sha256.lower())
    if not all(SHA256_PATTERN.fullmatch(value) for value in expected_hashes):
        raise StageAuditError("expected normalized SHA-256 values must be lowercase hex")
    if expected_pairs <= 0 or expected_length <= 0:
        raise StageAuditError("expected pair count and read length must be positive")

    normalized_paths = (normalized_r1, normalized_r2)
    initial_state = tuple(_stable_file_state(path) for path in normalized_paths)
    observed_hashes = tuple(_sha256(path) for path in normalized_paths)
    expected_sizes = (expected_r1_bytes, expected_r2_bytes)
    observed_sizes = tuple(state[0] for state in initial_state)
    if observed_hashes != expected_hashes or observed_sizes != expected_sizes:
        raise StageAuditError("normalized inputs differ from the pre-filter baseline")

    normalized = inspect_fastq_pair(
        normalized_r1,
        normalized_r2,
        expected_pairs=expected_pairs,
        expected_length=expected_length,
    )
    source = inspect_fastq_pair(
        source_r1,
        source_r2,
        expected_length=expected_length,
    )
    strict = inspect_fastq_pair(
        strict_r1,
        strict_r2,
        expected_length=expected_length,
    )
    final_state = tuple(_stable_file_state(path) for path in normalized_paths)
    if final_state != initial_state:
        raise StageAuditError("normalized inputs changed during the audit")

    return {
        "schema_version": 1,
        "checkpoint": "normalized_grch38_stage_audit",
        "captured_at_utc": datetime.now(UTC).isoformat(),
        "status": "PASS",
        "sample_id": sample_id,
        "normalized": {
            **normalized,
            "r1_bytes": observed_sizes[0],
            "r2_bytes": observed_sizes[1],
            "r1_sha256": observed_hashes[0],
            "r2_sha256": observed_hashes[1],
            "matches_prefilter_baseline": True,
            "stable_during_audit": True,
        },
        "grch38_filter": {
            "source": {
                **source,
                **_removal_metrics(expected_pairs, int(source["pairs"])),
            },
            "strict": {
                **strict,
                **_removal_metrics(expected_pairs, int(strict["pairs"])),
            },
        },
        "privacy": {
            "aggregate_evidence_only": True,
            "absolute_paths_recorded": False,
            "read_identifiers_recorded": False,
            "sequence_or_quality_recorded": False,
        },
    }
