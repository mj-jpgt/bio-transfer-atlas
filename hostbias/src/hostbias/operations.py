"""Safe production overlays and aggregate-only workflow progress."""

from __future__ import annotations

import hashlib
import os
import re
import shutil
import stat
import subprocess
import tempfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable

import yaml

from hostbias.config import load_and_validate
from hostbias.provenance import write_json_atomic


RUN_ID_PATTERN = re.compile(r"^[a-z0-9][a-z0-9._-]{2,80}$")
SHARED_FILESYSTEMS = {"nfs", "nfs4", "virtiofs"}
DEFAULT_SCHEDULER = {
    "cores": 30,
    "jobs": 8,
    "mem_mb": 204_800,
    "disk_mb": 2_500_000,
    "latency_wait_seconds": 120,
}
STAGE_TARGETS = {
    "fetch": "fetch_stage",
    "normalize": "normalize_stage",
    "filter": "filter_stage",
    "assemble": "assemble_stage",
    "downstream": "downstream_bridge",
    "endpoint": "gate_a_endpoint_aggregates",
}


class OperationError(ValueError):
    """Raised when a production operation would violate its safety contract."""


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _atomic_yaml(payload: dict[str, Any], output: Path) -> str:
    content = yaml.safe_dump(payload, sort_keys=False).encode("utf-8")
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, output)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise
    return hashlib.sha256(content).hexdigest()


def detect_mount_fstype(path: Path) -> str | None:
    """Return the Linux mount type for the longest mountpoint containing path."""

    mountinfo = Path("/proc/self/mountinfo")
    if not mountinfo.is_file():
        return None
    resolved = path.resolve()
    matches: list[tuple[int, str]] = []
    for line in mountinfo.read_text(encoding="utf-8").splitlines():
        left, separator, right = line.partition(" - ")
        if not separator:
            continue
        fields = left.split()
        right_fields = right.split()
        if len(fields) < 5 or not right_fields:
            continue
        mountpoint = Path(fields[4].replace("\\040", " "))
        try:
            resolved.relative_to(mountpoint)
        except ValueError:
            continue
        matches.append((len(mountpoint.parts), right_fields[0]))
    return max(matches)[1] if matches else None


def _safe_nfs_root(path: Path, *, require_shared_fs: bool) -> tuple[Path, str | None]:
    if not path.is_absolute():
        raise OperationError("NFS root must be an absolute path")
    if not path.exists() or not path.is_dir():
        raise OperationError("NFS root must already exist; refusing to create a mount path")
    resolved = path.resolve()
    if resolved == Path(resolved.anchor):
        raise OperationError("NFS root cannot be the filesystem root")
    filesystem = detect_mount_fstype(resolved)
    if require_shared_fs and filesystem not in SHARED_FILESYSTEMS:
        raise OperationError(
            f"NFS root filesystem must be one of {sorted(SHARED_FILESYSTEMS)}; "
            f"observed {filesystem or 'unknown'}"
        )
    return resolved, filesystem


def _under(path: Path, parent: Path, label: str) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(parent.resolve())
    except ValueError as error:
        raise OperationError(f"{label} must remain under {parent}") from error
    return resolved


def _nfs_reference(value: str, nfs_root: Path) -> str:
    path = Path(value)
    resolved = path.resolve() if path.is_absolute() else (nfs_root / path).resolve()
    _under(resolved, nfs_root, "reference path")
    return str(resolved)


def prepare_production_overlay(
    *,
    base_config: Path,
    output_config: Path,
    evidence_output: Path,
    nfs_root: Path,
    run_id: str,
    expected_samples: int = 40,
    require_shared_fs: bool = True,
) -> dict[str, Any]:
    """Write a production config whose sensitive work paths stay on shared NFS."""

    if not RUN_ID_PATTERN.fullmatch(run_id):
        raise OperationError("run_id must be 3-81 lowercase safe filename characters")
    inputs = load_and_validate(base_config)
    if len(inputs.samples) != expected_samples:
        raise OperationError(
            f"production config requires {expected_samples} samples; "
            f"observed {len(inputs.samples)}"
        )
    project_root = inputs.root.resolve()
    output_config = _under(output_config, project_root / "runtime", "output config")
    resolved_nfs, filesystem = _safe_nfs_root(
        nfs_root, require_shared_fs=require_shared_fs
    )

    scratch_root = _under(
        resolved_nfs / "scratch" / run_id,
        resolved_nfs,
        "scratch root",
    )
    work_dir = scratch_root / "work"
    work_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
    scratch_root.chmod(0o700)
    work_dir.chmod(0o700)
    scratch_mode = stat.S_IMODE(scratch_root.stat().st_mode)
    work_mode = stat.S_IMODE(work_dir.stat().st_mode)
    if os.name == "posix" and (scratch_mode != 0o700 or work_mode != 0o700):
        raise OperationError("scratch directories must have exact mode 0700")

    config = yaml.safe_load(base_config.read_text(encoding="utf-8"))
    config["paths"]["work_dir"] = str(work_dir)
    config["paths"]["results_dir"] = str(project_root / "results")
    config["references"]["grch38"]["fasta"] = _nfs_reference(
        config["references"]["grch38"]["fasta"], resolved_nfs
    )
    config["references"]["grch38"]["bowtie2_index"] = _nfs_reference(
        config["references"]["grch38"]["bowtie2_index"], resolved_nfs
    )
    for domain in ("human", "gtdb"):
        config["references"]["competitive"][domain]["minimap2_index"] = (
            _nfs_reference(
                config["references"]["competitive"][domain]["minimap2_index"],
                resolved_nfs,
            )
        )
    for name, resources in config["resources"].items():
        if resources["threads"] > DEFAULT_SCHEDULER["cores"]:
            raise OperationError(f"{name}: threads exceed production core cap")
        if resources["mem_mb"] > DEFAULT_SCHEDULER["mem_mb"]:
            raise OperationError(f"{name}: memory exceeds production memory cap")

    config_sha = _atomic_yaml(config, output_config)
    # Validate the generated overlay through the same workflow-facing contract.
    generated = load_and_validate(output_config)
    if len(generated.samples) != expected_samples:
        raise AssertionError("generated overlay changed sample count")

    evidence: dict[str, Any] = {
        "schema_version": 1,
        "status": "READY",
        "run_id": run_id,
        "sample_count": len(generated.samples),
        "filter_modes": list(config["filtering"]["modes"]),
        "base_config_sha256": sha256_file(base_config),
        "production_config_sha256": config_sha,
        "sample_manifest_sha256": sha256_file(generated.manifest_path),
        "storage": {
            "shared_filesystem_type": filesystem,
            "scratch_layout": f"scratch/{run_id}/work",
            "scratch_mode": "0700",
            "sequence_data_in_git_results": False,
        },
        "scheduler": dict(DEFAULT_SCHEDULER),
        "privacy": {
            "contains_absolute_paths": False,
            "contains_fastq_urls": False,
            "contains_fastq_checksums": False,
            "contains_sequence_data": False,
        },
    }
    write_json_atomic(evidence, evidence_output)
    return evidence


def _complete(paths: tuple[Path, ...]) -> bool:
    return all(path.exists() and path.is_file() for path in paths)


def production_status(config_path: Path) -> dict[str, Any]:
    """Summarize stage completion without publishing paths or sample-level state."""

    inputs = load_and_validate(config_path)
    samples = [row["sample_id"] for row in inputs.samples if row["status"] == "primary"]
    modes = tuple(inputs.config["filtering"]["modes"])
    work = Path(inputs.config["paths"]["work_dir"])
    results = Path(inputs.config["paths"]["results_dir"])

    units: dict[str, list[tuple[Path, ...]]] = {
        "fetch": [
            (work / "raw" / f"{sample}_R1.fastq.gz", work / "raw" / f"{sample}_R2.fastq.gz")
            for sample in samples
        ],
        "normalize": [
            (
                work / "normalized" / f"{sample}_R1.fastq.gz",
                work / "normalized" / f"{sample}_R2.fastq.gz",
                work / "normalized" / f"{sample}.fastp.json",
            )
            for sample in samples
        ],
        "filter": [
            (
                work / "filtered" / sample / "source_R1.fastq.gz",
                work / "filtered" / sample / "source_R2.fastq.gz",
                work / "filtered" / sample / "strict_R1.fastq.gz",
                work / "filtered" / sample / "strict_R2.fastq.gz",
                work / "filtered" / sample / "bowtie2.metrics.txt",
            )
            for sample in samples
        ],
        "assemble": [
            (
                work / "assembly" / sample / mode / "final.contigs.fa",
                work / "assembly" / sample / mode / "megahit.log",
            )
            for sample in samples
            for mode in modes
        ],
        "downstream": [
            (
                results / "aggregate" / "assembly_qc" / f"{sample}.{mode}.json",
                results
                / "aggregate"
                / "alignment_manifests"
                / f"{sample}.{mode}.json",
            )
            for sample in samples
            for mode in modes
        ],
        "endpoint": [
            (
                results
                / "aggregate"
                / "sample_endpoints"
                / f"{sample}.{mode}.json",
            )
            for sample in samples
            for mode in modes
        ],
    }
    stages: dict[str, dict[str, Any]] = {}
    for stage, stage_units in units.items():
        completed = sum(_complete(unit) for unit in stage_units)
        expected = len(stage_units)
        state = "complete" if completed == expected else ("not_started" if completed == 0 else "partial")
        stages[stage] = {
            "state": state,
            "completed_units": completed,
            "expected_units": expected,
            "percent": round(100 * completed / expected, 2) if expected else 100.0,
        }

    scratch_root = work.parent
    mode = stat.S_IMODE(scratch_root.stat().st_mode) if scratch_root.exists() else None
    disk = shutil.disk_usage(work if work.exists() else scratch_root.parent)
    return {
        "schema_version": 1,
        "observed_at_utc": datetime.now(UTC).isoformat(),
        "experiment_id": inputs.config["experiment"]["id"],
        "sample_count": len(samples),
        "filter_modes": list(modes),
        "config_sha256": sha256_file(inputs.config_path),
        "stages": stages,
        "storage": {
            "scratch_exists": scratch_root.is_dir(),
            "scratch_mode_ok": mode == 0o700 if os.name == "posix" else None,
            "free_bytes": disk.free,
            "total_bytes": disk.total,
        },
        "privacy": {
            "contains_absolute_paths": False,
            "contains_sample_accessions": False,
            "contains_sequence_data": False,
        },
    }


def _git_value(root: Path, *arguments: str) -> str | None:
    completed = subprocess.run(
        ["git", *arguments],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def build_snakemake_command(
    *,
    config_path: Path,
    stage: str,
    cores: int = 30,
    jobs: int = 8,
    mem_mb: int = 204_800,
    disk_mb: int = 2_500_000,
    latency_wait_seconds: int = 120,
    snakemake_executable: str = "snakemake",
    dry_run: bool = False,
) -> list[str]:
    if stage not in STAGE_TARGETS:
        raise OperationError(f"stage must be one of {sorted(STAGE_TARGETS)}")
    if not 1 <= cores <= DEFAULT_SCHEDULER["cores"]:
        raise OperationError("cores must be between 1 and the 30-core VM cap")
    if not 1 <= jobs <= DEFAULT_SCHEDULER["jobs"]:
        raise OperationError("jobs must be between 1 and the production cap of 8")
    if not 1 <= mem_mb <= DEFAULT_SCHEDULER["mem_mb"]:
        raise OperationError("mem_mb exceeds the production memory cap")
    if not 1 <= disk_mb <= DEFAULT_SCHEDULER["disk_mb"]:
        raise OperationError("disk_mb exceeds the production disk cap")
    command = [
        snakemake_executable,
        STAGE_TARGETS[stage],
        "--snakefile",
        "workflow/Snakefile",
        "--configfile",
        str(config_path.resolve()),
        "--cores",
        str(cores),
        "--jobs",
        str(jobs),
        "--resources",
        f"mem_mb={mem_mb}",
        f"disk_mb={disk_mb}",
        "--rerun-incomplete",
        "--keep-going",
        "--latency-wait",
        str(latency_wait_seconds),
        "--use-conda",
        "--conda-frontend",
        "conda",
        "--printshellcmds",
        "--rerun-triggers",
        "mtime",
        "input",
        "params",
        "code",
        "software-env",
    ]
    if dry_run:
        command.append("--dry-run")
    return command


Runner = Callable[[list[str], Path], int]


def _default_runner(command: list[str], cwd: Path) -> int:
    return subprocess.run(command, cwd=cwd, check=False).returncode


def launch_production(
    *,
    config_path: Path,
    stage: str,
    evidence_output: Path,
    cores: int = 30,
    jobs: int = 8,
    mem_mb: int = 204_800,
    disk_mb: int = 2_500_000,
    latency_wait_seconds: int = 120,
    snakemake_executable: str = "snakemake",
    dry_run: bool = False,
    require_clean_git: bool = True,
    require_shared_fs: bool = True,
    runner: Runner | None = None,
) -> dict[str, Any]:
    """Run one production stage in the foreground and record restart evidence."""

    inputs = load_and_validate(config_path)
    if len(inputs.samples) != 40:
        raise OperationError(
            f"production launch requires 40 samples; observed {len(inputs.samples)}"
        )
    work_dir = Path(inputs.config["paths"]["work_dir"])
    results_dir = Path(inputs.config["paths"]["results_dir"])
    if require_shared_fs:
        if not work_dir.is_absolute() or not work_dir.is_dir():
            raise OperationError("production work directory must already exist on NFS")
        filesystem = detect_mount_fstype(work_dir)
        if filesystem not in SHARED_FILESYSTEMS:
            raise OperationError(
                "production work directory is not on an approved shared filesystem"
            )
        if stat.S_IMODE(work_dir.parent.stat().st_mode) != 0o700:
            raise OperationError("production scratch root must have exact mode 0700")
        if results_dir.resolve() != (inputs.root / "results").resolve():
            raise OperationError(
                "production aggregate results must stay in the project results tree"
            )
        references = [
            inputs.config["references"]["grch38"]["fasta"],
            inputs.config["references"]["grch38"]["bowtie2_index"],
            inputs.config["references"]["competitive"]["human"]["minimap2_index"],
            inputs.config["references"]["competitive"]["gtdb"]["minimap2_index"],
        ]
        if any(
            not Path(reference).is_absolute()
            or detect_mount_fstype(Path(reference).parent) not in SHARED_FILESYSTEMS
            for reference in references
        ):
            raise OperationError(
                "all production references must resolve to shared storage"
            )
    tracked_status = _git_value(
        inputs.root, "status", "--porcelain", "--untracked-files=no"
    )
    if require_clean_git and tracked_status:
        raise OperationError("production launch requires a clean tracked worktree")
    commit = _git_value(inputs.root, "rev-parse", "HEAD")
    command = build_snakemake_command(
        config_path=config_path,
        stage=stage,
        cores=cores,
        jobs=jobs,
        mem_mb=mem_mb,
        disk_mb=disk_mb,
        latency_wait_seconds=latency_wait_seconds,
        snakemake_executable=snakemake_executable,
        dry_run=dry_run,
    )
    started = datetime.now(UTC).isoformat()
    evidence: dict[str, Any] = {
        "schema_version": 1,
        "status": "RUNNING",
        "stage": stage,
        "target": STAGE_TARGETS[stage],
        "started_at_utc": started,
        "finished_at_utc": None,
        "exit_code": None,
        "dry_run": dry_run,
        "git_commit": commit,
        "config_sha256": sha256_file(inputs.config_path),
        "scheduler": {
            "cores": cores,
            "jobs": jobs,
            "mem_mb": mem_mb,
            "disk_mb": disk_mb,
            "latency_wait_seconds": latency_wait_seconds,
        },
        "resume": {
            "rerun_incomplete": True,
            "keep_going": True,
            "force_outputs": False,
            "delete_outputs": False,
            "same_command_resumes": True,
        },
        "privacy": {
            "contains_absolute_paths": False,
            "contains_sample_accessions": False,
            "contains_sequence_data": False,
        },
    }
    write_json_atomic(evidence, evidence_output)
    exit_code = (runner or _default_runner)(command, inputs.root)
    evidence["exit_code"] = exit_code
    evidence["finished_at_utc"] = datetime.now(UTC).isoformat()
    evidence["status"] = (
        "DRY_RUN_COMPLETE"
        if dry_run and exit_code == 0
        else ("COMPLETE" if exit_code == 0 else "FAILED")
    )
    evidence["progress"] = production_status(config_path)["stages"]
    write_json_atomic(evidence, evidence_output)
    return evidence
