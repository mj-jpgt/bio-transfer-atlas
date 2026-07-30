"""Generate privacy-conscious, machine-readable provenance."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import tempfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from hostbias.config import ValidatedInputs


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git(root: Path, *args: str) -> str | None:
    completed = subprocess.run(
        ["git", *args],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def build_provenance(inputs: ValidatedInputs) -> dict[str, Any]:
    """Return aggregate-safe provenance without environment or absolute data paths."""

    commit = _git(inputs.root, "rev-parse", "HEAD")
    status = _git(inputs.root, "status", "--porcelain", "--untracked-files=no")
    return {
        "schema_version": 1,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "experiment_id": inputs.config["experiment"]["id"],
        "seed": inputs.config["experiment"]["seed"],
        "git": {
            "commit": commit,
            "dirty_tracked_files": bool(status),
        },
        "runtime": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "platform": platform.platform(),
        },
        "inputs": {
            "config": {
                "path": inputs.config_path.relative_to(inputs.root).as_posix(),
                "sha256": sha256_file(inputs.config_path),
            },
            "sample_manifest": {
                "path": inputs.manifest_path.relative_to(inputs.root).as_posix(),
                "sha256": sha256_file(inputs.manifest_path),
                "sample_count": len(inputs.samples),
            },
        },
        "privacy": {
            "environment_variables_recorded": False,
            "absolute_data_paths_recorded": False,
            "read_level_data_recorded": False,
        },
    }


def write_json_atomic(payload: dict[str, Any], output: Path) -> None:
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, output)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise
