"""Safe production overlays and aggregate-only workflow progress."""

from __future__ import annotations

import hashlib
import os
import re
import stat
import tempfile
from pathlib import Path
from typing import Any

import yaml

from hostbias.config import ValidationError, load_and_validate
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
