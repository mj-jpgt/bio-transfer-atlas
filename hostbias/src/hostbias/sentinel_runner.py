"""Restartable, sequence-discarding execution of the six Gate A sentinels."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import time
from contextlib import ExitStack
from datetime import UTC, datetime
from pathlib import Path
from typing import BinaryIO, Mapping, Protocol, Sequence, TextIO

import yaml

from .data_manifest import (
    MANIFEST_FIELDS,
    ManifestError,
    canonical_tsv,
    read_tsv,
    sha256_bytes,
    validate_manifest,
)
from .provenance import sha256_file, write_json_atomic
from .sentinel import evaluate_sentinels

ACCESSION_PATTERN = re.compile(r"^[A-Z0-9]+$")
SEQUENCE_PATTERN = re.compile(r"(?i)\b[acgtn]{20,}\b")
SENTINEL_REPORT = "sentinel_eligibility.json"


class SentinelExecutionError(RuntimeError):
    """Operational error carrying a sequence-safe stage diagnostic."""

    def __init__(
        self,
        stage: str,
        message: str,
        *,
        return_code: int | None = None,
        stderr: bytes = b"",
    ) -> None:
        super().__init__(message)
        self.stage = stage
        self.return_code = return_code
        self.stderr = stderr


class Executor(Protocol):
    def tool_versions(self) -> dict[str, str]: ...

    def fetch(
        self, accession: str, scratch: Path, spots: int, threads: int
    ) -> tuple[Path, Path]: ...

    def align(
        self,
        r1: Path,
        r2: Path,
        grch38_index: Path,
        threads: int,
        scratch: Path,
    ) -> int: ...


def _tool_version(executable: str) -> str:
    try:
        completed = subprocess.run(
            [executable, "--version"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise SentinelExecutionError("preflight", f"cannot execute {executable}") from exc
    output = completed.stdout.decode("utf-8", errors="replace")
    first_line = next((line.strip() for line in output.splitlines() if line.strip()), "")
    if completed.returncode != 0 or not first_line:
        raise SentinelExecutionError(
            "preflight",
            f"{executable} --version failed",
            return_code=completed.returncode,
            stderr=completed.stdout,
        )
    return first_line[:300]


def consume_mapped_pair_sam(stream: BinaryIO) -> int:
    """Count query-name groups with a mapped primary alignment from ordered SAM."""
    mapped_pairs = 0
    previous_qname: bytes | None = None
    previous_counted = False
    for line in stream:
        if line.startswith(b"@"):
            continue
        fields = line.split(b"\t", 2)
        if len(fields) < 2:
            raise SentinelExecutionError("alignment", "bowtie2 emitted malformed SAM")
        qname = fields[0]
        try:
            flag = int(fields[1])
        except ValueError as exc:
            raise SentinelExecutionError(
                "alignment", "bowtie2 emitted a non-integer SAM flag"
            ) from exc
        if qname != previous_qname:
            previous_qname = qname
            previous_counted = False
        primary_mapped = not flag & (0x4 | 0x100 | 0x800)
        if primary_mapped and not previous_counted:
            mapped_pairs += 1
            previous_counted = True
    return mapped_pairs


class LocalExecutor:
    """Invoke SRA Toolkit and Bowtie2 without materializing alignment output."""

    def __init__(
        self,
        *,
        fastq_dump: str = "fastq-dump",
        bowtie2: str = "bowtie2",
    ) -> None:
        self.fastq_dump = fastq_dump
        self.bowtie2 = bowtie2

    def tool_versions(self) -> dict[str, str]:
        return {
            "fastq_dump": _tool_version(self.fastq_dump),
            "bowtie2": _tool_version(self.bowtie2),
        }

    def fetch(
        self, accession: str, scratch: Path, spots: int, threads: int
    ) -> tuple[Path, Path]:
        command = [
            self.fastq_dump,
            "--split-files",
            "-O",
            str(scratch),
            "-N",
            "1",
            "-X",
            str(spots),
            accession,
        ]
        try:
            completed = subprocess.run(
                command,
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
            )
        except OSError as exc:
            raise SentinelExecutionError(
                "fetch", "fastq-dump could not start"
            ) from exc
        if completed.returncode:
            raise SentinelExecutionError(
                "fetch",
                "fastq-dump failed",
                return_code=completed.returncode,
                stderr=completed.stderr,
            )
        r1 = scratch / f"{accession}_1.fastq"
        r2 = scratch / f"{accession}_2.fastq"
        observed_fastqs = sorted(scratch.glob(f"{accession}*.fastq"))
        if not r1.is_file() or not r2.is_file() or observed_fastqs != [r1, r2]:
            raise SentinelExecutionError(
                "pair_validation",
                "fastq-dump did not produce exactly two mate FASTQs",
            )
        return r1, r2

    def align(
        self,
        r1: Path,
        r2: Path,
        grch38_index: Path,
        threads: int,
        scratch: Path,
    ) -> int:
        stderr_path = scratch / "bowtie2.stderr"
        command = [
            self.bowtie2,
            "--very-sensitive-local",
            "--reorder",
            "--no-unal",
            "-x",
            str(grch38_index),
            "-1",
            str(r1),
            "-2",
            str(r2),
            "-p",
            str(threads),
        ]
        try:
            with stderr_path.open("wb") as stderr_handle:
                process = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=stderr_handle,
                )
                assert process.stdout is not None
                try:
                    mapped_pairs = consume_mapped_pair_sam(process.stdout)
                except BaseException:
                    process.kill()
                    process.wait()
                    raise
                return_code = process.wait()
        except OSError as exc:
            raise SentinelExecutionError(
                "alignment", "bowtie2 could not start"
            ) from exc
        stderr = stderr_path.read_bytes() if stderr_path.exists() else b""
        if return_code:
            raise SentinelExecutionError(
                "alignment",
                "bowtie2 failed",
                return_code=return_code,
                stderr=stderr,
            )
        return mapped_pairs


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii")
    return path.open("r", encoding="ascii")


def _canonical_fastq_name(header: str, record_number: int) -> str:
    if not header.startswith("@"):
        raise SentinelExecutionError(
            "pair_validation", f"record {record_number} has an invalid FASTQ header"
        )
    name = header[1:].split()[0]
    if name.endswith("/1") or name.endswith("/2"):
        return name[:-2]
    return name


def _read_fastq_record(
    handle: TextIO, record_number: int
) -> tuple[str, int] | None:
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline().rstrip("\r\n")
    separator = handle.readline()
    quality = handle.readline().rstrip("\r\n")
    if not separator.startswith("+") or len(sequence) != len(quality):
        raise SentinelExecutionError(
            "pair_validation", f"record {record_number} has invalid FASTQ shape"
        )
    return _canonical_fastq_name(header, record_number), len(sequence)


def validate_fastq_pair(r1: Path, r2: Path) -> tuple[int, dict[str, int]]:
    """Validate mate synchronization and return only aggregate shape counts."""
    count = 0
    minimum_length: int | None = None
    maximum_length = 0
    total_bases = 0
    with ExitStack() as stack:
        first = stack.enter_context(_open_text(r1))
        second = stack.enter_context(_open_text(r2))
        while True:
            record_number = count + 1
            first_record = _read_fastq_record(first, record_number)
            second_record = _read_fastq_record(second, record_number)
            if first_record is None and second_record is None:
                break
            if first_record is None or second_record is None:
                raise SentinelExecutionError(
                    "pair_validation", "mate FASTQs contain different record counts"
                )
            if first_record[0] != second_record[0]:
                raise SentinelExecutionError(
                    "pair_validation", f"mate names differ at record {record_number}"
                )
            lengths = (first_record[1], second_record[1])
            minimum_length = min(
                lengths if minimum_length is None else (*lengths, minimum_length)
            )
            maximum_length = max(maximum_length, *lengths)
            total_bases += sum(lengths)
            count += 1
    if count == 0:
        raise SentinelExecutionError("pair_validation", "mate FASTQs are empty")
    assert minimum_length is not None
    return count, {
        "minimum_read_length": minimum_length,
        "maximum_read_length": maximum_length,
        "total_bases": total_bases,
    }


def fingerprint_bowtie2_index(prefix: Path) -> dict[str, object]:
    candidates = sorted(
        [
            *prefix.parent.glob(f"{prefix.name}*.bt2"),
            *prefix.parent.glob(f"{prefix.name}*.bt2l"),
        ],
        key=lambda path: path.name,
    )
    if not candidates:
        raise SentinelExecutionError(
            "preflight", "no Bowtie2 index files found for configured GRCh38 prefix"
        )
    files = [{"name": path.name, "sha256": sha256_file(path)} for path in candidates]
    combined = hashlib.sha256()
    for item in files:
        combined.update(item["name"].encode("utf-8"))
        combined.update(b"\0")
        combined.update(item["sha256"].encode("ascii"))
        combined.update(b"\n")
    return {"sha256": combined.hexdigest(), "files": files}


def _safe_diagnostic_text(payload: bytes) -> tuple[str, int, str]:
    digest = hashlib.sha256(payload).hexdigest()
    decoded = payload.decode("utf-8", errors="replace")
    safe_lines = []
    for line in decoded.splitlines()[-30:]:
        line = "".join(character for character in line if character.isprintable())
        line = SEQUENCE_PATTERN.sub("<redacted-sequence>", line)
        safe_lines.append(line[:500])
    return "\n".join(safe_lines), len(payload), digest


def _selected_sentinels(
    manifest: Sequence[Mapping[str, str]], runs_per_arm: int
) -> list[Mapping[str, str]]:
    selected = [
        row
        for row in manifest
        if row["role"] == "primary" and int(row["rank"]) <= runs_per_arm
    ]
    arms = sorted({row["arm"] for row in manifest})
    for arm in arms:
        arm_rows = [row for row in selected if row["arm"] == arm]
        if [int(row["rank"]) for row in arm_rows] != list(range(1, runs_per_arm + 1)):
            raise ManifestError(f"{arm}: frozen sentinels must be primary ranks 1..{runs_per_arm}")
    return selected


def _checkpoint_matches(
    checkpoint: Mapping[str, object], expected: Mapping[str, object]
) -> bool:
    return checkpoint.get("status") == "complete" and all(
        checkpoint.get(key) == value for key, value in expected.items()
    )


def _load_json(path: Path) -> dict[str, object] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def _remove_run_scratch(run_scratch: Path, scratch_root: Path) -> None:
    resolved = run_scratch.resolve()
    root = scratch_root.resolve()
    if resolved.parent != root or not ACCESSION_PATTERN.fullmatch(resolved.name):
        raise SentinelExecutionError("cleanup", "refusing unsafe scratch cleanup target")
    if resolved.exists():
        shutil.rmtree(resolved)


def _run_one(
    sample: Mapping[str, str],
    *,
    spots: int,
    threads: int,
    manifest_sha256: str,
    reference: Mapping[str, object],
    tools: Mapping[str, str],
    grch38_index: Path,
    scratch_root: Path,
    output_dir: Path,
    executor: Executor,
) -> dict[str, object]:
    accession = sample["run_accession"]
    if not ACCESSION_PATTERN.fullmatch(accession):
        raise ManifestError(f"unsafe run accession: {accession!r}")
    expected = {
        "run_accession": accession,
        "arm": sample["arm"],
        "rank": int(sample["rank"]),
        "requested_spots": spots,
        "manifest_sha256": manifest_sha256,
        "grch38_index_sha256": reference["sha256"],
        "tool_versions": dict(tools),
    }
    checkpoint_path = output_dir / "runs" / f"{accession}.json"
    checkpoint = _load_json(checkpoint_path)
    if checkpoint and _checkpoint_matches(checkpoint, expected):
        return checkpoint

    run_scratch = scratch_root / accession
    _remove_run_scratch(run_scratch, scratch_root)
    run_scratch.mkdir(parents=True)
    diagnostic_path = output_dir / "diagnostics" / f"{accession}.json"
    start_total = time.monotonic()
    try:
        start = time.monotonic()
        r1, r2 = executor.fetch(accession, run_scratch, spots, threads)
        fetch_seconds = time.monotonic() - start

        start = time.monotonic()
        observed_pairs, shape = validate_fastq_pair(r1, r2)
        checksum_before = {"r1": sha256_file(r1), "r2": sha256_file(r2)}
        validation_seconds = time.monotonic() - start

        start = time.monotonic()
        mapped_pairs = executor.align(
            r1, r2, grch38_index, threads, run_scratch
        )
        alignment_seconds = time.monotonic() - start
        if mapped_pairs < 0 or mapped_pairs > observed_pairs:
            raise SentinelExecutionError(
                "alignment", "mapped-pair count is outside the observed pair count"
            )
        checksum_after = {"r1": sha256_file(r1), "r2": sha256_file(r2)}
        checksums_stable = checksum_before == checksum_after
        if not checksums_stable:
            raise SentinelExecutionError(
                "checksum", "input FASTQ checksum changed during alignment"
            )
        result: dict[str, object] = {
            "schema_version": 1,
            "status": "complete",
            **expected,
            "observed_pairs": observed_pairs,
            "grch38_mapped_pairs": mapped_pairs,
            "mapped_pairs_per_million": mapped_pairs / observed_pairs * 1_000_000,
            "fastq_shape": shape,
            "input_sha256": checksum_before,
            "checksums_stable": True,
            "metadata_ok": True,
            "checksum_ok": True,
            "timing_seconds": {
                "fetch": round(fetch_seconds, 6),
                "pair_validation_and_checksum": round(validation_seconds, 6),
                "alignment": round(alignment_seconds, 6),
                "total": round(time.monotonic() - start_total, 6),
            },
        }
        _remove_run_scratch(run_scratch, scratch_root)
        write_json_atomic(result, checkpoint_path)
        diagnostic_path.unlink(missing_ok=True)
        return result
    except Exception as error:
        if isinstance(error, SentinelExecutionError):
            stage = error.stage
            return_code = error.return_code
            safe_tail, stderr_bytes, stderr_sha256 = _safe_diagnostic_text(error.stderr)
            message = str(error)
        else:
            stage = "unexpected"
            return_code = None
            safe_tail, stderr_bytes, stderr_sha256 = "", 0, hashlib.sha256(b"").hexdigest()
            message = type(error).__name__
        scratch_removed = False
        cleanup_error: SentinelExecutionError | None = None
        try:
            _remove_run_scratch(run_scratch, scratch_root)
            scratch_removed = True
        except SentinelExecutionError as cleanup_failure:
            cleanup_error = cleanup_failure
        diagnostic = {
            "schema_version": 1,
            "status": "failed",
            "run_accession": accession,
            "arm": sample["arm"],
            "rank": int(sample["rank"]),
            "stage": stage,
            "message": message[:500],
            "return_code": return_code,
            "stderr_bytes": stderr_bytes,
            "stderr_sha256": stderr_sha256,
            "sanitized_stderr_tail": safe_tail,
            "sequence_data_recorded": False,
            "scratch_removed": scratch_removed,
        }
        if cleanup_error is not None:
            diagnostic["cleanup_error"] = str(cleanup_error)
        write_json_atomic(diagnostic, diagnostic_path)
        if cleanup_error is not None:
            raise cleanup_error from error
        raise


def run_sentinel_panel(
    *,
    manifest_path: Path,
    thresholds_path: Path,
    grch38_index: Path,
    scratch_root: Path,
    output_dir: Path,
    threads: int,
    executor: Executor | None = None,
) -> dict[str, object]:
    """Run or resume all frozen sentinels and write an aggregate-only verdict."""
    manifest = read_tsv(manifest_path)
    validate_manifest(manifest)
    thresholds = yaml.safe_load(thresholds_path.read_text(encoding="utf-8"))
    sentinel_thresholds = thresholds["sentinel"]
    runs_per_arm = int(sentinel_thresholds["runs_per_arm"])
    spots = int(sentinel_thresholds["streamed_spots_per_run"])
    minimum_rate = float(
        sentinel_thresholds["minimum_grch38_mapped_pairs_per_million"]
    )
    minimum_passing = int(sentinel_thresholds["minimum_passing_runs_per_arm"])
    samples = _selected_sentinels(manifest, runs_per_arm)
    manifest_sha256 = sha256_bytes(canonical_tsv(manifest, MANIFEST_FIELDS))

    executor = executor or LocalExecutor()
    tools = executor.tool_versions()
    reference = fingerprint_bowtie2_index(grch38_index)
    output_dir.mkdir(parents=True, exist_ok=True)
    scratch_root.mkdir(parents=True, exist_ok=True)
    completed: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    for sample in samples:
        try:
            completed.append(
                _run_one(
                    sample,
                    spots=spots,
                    threads=threads,
                    manifest_sha256=manifest_sha256,
                    reference=reference,
                    tools=tools,
                    grch38_index=grch38_index,
                    scratch_root=scratch_root,
                    output_dir=output_dir,
                    executor=executor,
                )
            )
        except Exception:
            diagnostic_path = (
                output_dir / "diagnostics" / f"{sample['run_accession']}.json"
            )
            diagnostic = _load_json(diagnostic_path) or {}
            failures.append(
                {
                    "run_accession": sample["run_accession"],
                    "stage": diagnostic.get("stage", "unexpected"),
                    "diagnostic": f"diagnostics/{sample['run_accession']}.json",
                }
            )

    report: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "status": "complete" if not failures else "operational_failure",
        "eligible": False,
        "manifest_sha256": manifest_sha256,
        "grch38_index": reference,
        "tool_versions": tools,
        "completed_runs": len(completed),
        "required_runs": len(samples),
        "failures": failures,
        "runs": completed,
        "privacy": {
            "sequence_data_recorded": False,
            "sam_or_bam_materialized": False,
            "scratch_cleaned_after_each_run": True,
        },
    }
    if not failures:
        metrics = [
            {
                "run_accession": str(run["run_accession"]),
                "arm": str(run["arm"]),
                "streamed_spots": str(run["observed_pairs"]),
                "grch38_mapped_pairs": str(run["grch38_mapped_pairs"]),
                "metadata_ok": str(run["metadata_ok"]).lower(),
                "checksum_ok": str(run["checksum_ok"]).lower(),
            }
            for run in completed
        ]
        decision = evaluate_sentinels(
            manifest,
            metrics,
            sentinels_per_arm=runs_per_arm,
            minimum_mapped_pairs_per_million=minimum_rate,
            minimum_passing_runs=minimum_passing,
        )
        report["eligible"] = decision["eligible"]
        report["decision"] = decision
    write_json_atomic(report, output_dir / SENTINEL_REPORT)
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--thresholds", required=True, type=Path)
    parser.add_argument("--grch38-index", required=True, type=Path)
    parser.add_argument("--scratch-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", default=16, type=int)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    report = run_sentinel_panel(
        manifest_path=args.manifest,
        thresholds_path=args.thresholds,
        grch38_index=args.grch38_index,
        scratch_root=args.scratch_root,
        output_dir=args.output_dir,
        threads=args.threads,
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "eligible": report["eligible"],
                "report": str(args.output_dir / SENTINEL_REPORT),
            }
        )
    )
    return 0 if report["status"] == "complete" else 2


if __name__ == "__main__":
    raise SystemExit(main())
