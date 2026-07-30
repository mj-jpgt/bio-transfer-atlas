"""Build executable workflow inputs from frozen selection and ENA snapshots."""

from __future__ import annotations

import hashlib
import json
import re
import urllib.parse
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

import yaml

from hostbias.data_manifest import (
    MANIFEST_FIELDS,
    ManifestError,
    RUNTIME_ENA_FIELDS,
    canonical_tsv,
    read_tsv,
    sha256_bytes,
    validate_manifest,
)
from hostbias.provenance import write_json_atomic


RUNTIME_MANIFEST_FIELDS = (
    "sample_id",
    "cohort",
    "accession",
    "layout",
    "fastq1_url",
    "fastq1_md5",
    "fastq1_bytes",
    "fastq2_url",
    "fastq2_md5",
    "fastq2_bytes",
    "status",
    "rank",
)
SCOPES = {"sentinel", "primary"}
MD5_PATTERN = re.compile(r"^[a-fA-F0-9]{32}$")
ALLOWED_REPLACEMENT_REASONS = {
    "checksum_failure",
    "wrong_library_layout",
    "corrupt_or_unsynchronized_pairs",
    "insufficient_clean_pairs",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _split_pair(value: str, *, run: str, field: str) -> tuple[str, str]:
    parts = tuple(part.strip() for part in value.split(";") if part.strip())
    if len(parts) != 2:
        raise ManifestError(f"{run}: {field} must contain exactly two values")
    return parts


def _url(value: str, *, run: str) -> str:
    candidate = value if "://" in value else f"ftp://{value}"
    parsed = urllib.parse.urlparse(candidate)
    if parsed.scheme not in {"ftp", "https"} or not parsed.netloc:
        raise ManifestError(f"{run}: invalid public FASTQ URL")
    if parsed.username or parsed.password or parsed.query or parsed.fragment:
        raise ManifestError(f"{run}: FASTQ URL must not contain credentials or tokens")
    return candidate


def _positive_size(value: str, *, run: str) -> str:
    try:
        size = int(value)
    except ValueError as error:
        raise ManifestError(f"{run}: invalid fastq_bytes value") from error
    if size <= 0:
        raise ManifestError(f"{run}: fastq_bytes values must be positive")
    return str(size)


def _replacement_rows(
    frozen_rows: Sequence[Mapping[str, str]],
    replacements: Sequence[Mapping[str, object]],
) -> dict[str, Mapping[str, str]]:
    """Validate technical substitutions without changing the frozen ranking."""

    by_accession = {row["run_accession"]: row for row in frozen_rows}
    selected: dict[str, Mapping[str, str]] = {}
    used_replacements: set[str] = set()
    replaced_accessions: set[str] = set()
    for replacement in replacements:
        required = {"arm", "replaced_accession", "replacement_accession", "reason"}
        missing = required.difference(replacement)
        if missing:
            raise ManifestError(f"replacement missing fields {sorted(missing)}")
        arm = str(replacement["arm"])
        replaced = str(replacement["replaced_accession"])
        candidate = str(replacement["replacement_accession"])
        reason = str(replacement["reason"])
        if reason not in ALLOWED_REPLACEMENT_REASONS:
            raise ManifestError(f"{replaced}: replacement reason is not technical")
        source = by_accession.get(replaced)
        reserve = by_accession.get(candidate)
        if source is None or reserve is None:
            raise ManifestError("replacement accession absent from frozen manifest")
        if source["arm"] != arm or reserve["arm"] != arm:
            raise ManifestError(f"{replaced}: replacement must remain within arm")
        if source["role"] != "primary" or reserve["role"] != "reserve":
            raise ManifestError(f"{replaced}: replacement must be primary-to-reserve")
        if replaced in replaced_accessions or candidate in used_replacements:
            raise ManifestError("replacement ledger contains duplicate accessions")
        earlier_available = [
            row
            for row in frozen_rows
            if row["arm"] == arm
            and row["role"] == "reserve"
            and row["run_accession"] not in used_replacements
        ]
        first_available = min(earlier_available, key=lambda row: int(row["rank"]))
        if candidate != first_available["run_accession"]:
            raise ManifestError(
                f"{replaced}: replacement must use the next reserve in frozen rank order"
            )
        selected[replaced] = reserve
        used_replacements.add(candidate)
        replaced_accessions.add(replaced)
    return selected


def _selected_rows(
    frozen_rows: Sequence[Mapping[str, str]],
    scope: str,
    replacements: Sequence[Mapping[str, object]] = (),
) -> list[Mapping[str, str]]:
    if scope not in SCOPES:
        raise ManifestError(f"scope must be one of {sorted(SCOPES)}")
    validate_manifest(frozen_rows)
    if scope == "sentinel":
        selected = [
            row
            for row in frozen_rows
            if row["role"] == "primary" and int(row["rank"]) <= 3
        ]
    else:
        selected = [row for row in frozen_rows if row["role"] == "primary"]
    expected_per_arm = 3 if scope == "sentinel" else 20
    counts = Counter(row["arm"] for row in selected)
    if any(value != expected_per_arm for value in counts.values()) or len(counts) != 2:
        raise ManifestError(
            f"{scope}: expected two arms with {expected_per_arm} runs each; "
            f"observed {dict(counts)}"
        )
    if scope == "primary" and replacements:
        substitution = _replacement_rows(frozen_rows, replacements)
        selected = [substitution.get(row["run_accession"], row) for row in selected]
    return selected


def build_runtime_rows(
    frozen_rows: Sequence[Mapping[str, str]],
    snapshots: Mapping[str, Sequence[Mapping[str, str]]],
    *,
    scope: str,
    replacements: Sequence[Mapping[str, object]] = (),
) -> list[dict[str, str]]:
    """Join a frozen scope to checksum/size-pinned ENA FASTQ metadata."""

    selected = _selected_rows(frozen_rows, scope, replacements)
    frozen_arms = {row["arm"] for row in frozen_rows}
    if set(snapshots) != frozen_arms:
        raise ManifestError(
            f"snapshot arms must exactly match frozen arms {sorted(frozen_arms)}"
        )

    indexes: dict[str, dict[str, Mapping[str, str]]] = {}
    all_snapshot_runs: list[str] = []
    for arm in sorted(snapshots):
        rows = snapshots[arm]
        if not rows:
            raise ManifestError(f"{arm}: ENA snapshot is empty")
        missing = set(RUNTIME_ENA_FIELDS).difference(rows[0])
        if missing:
            raise ManifestError(f"{arm}: ENA snapshot missing fields {sorted(missing)}")
        runs = [row["run_accession"] for row in rows]
        if len(runs) != len(set(runs)):
            raise ManifestError(f"{arm}: ENA snapshot has duplicate run accessions")
        if runs != sorted(runs):
            raise ManifestError(f"{arm}: ENA snapshot must be canonical accession order")
        indexes[arm] = {row["run_accession"]: row for row in rows}
        all_snapshot_runs.extend(runs)
    if len(all_snapshot_runs) != len(set(all_snapshot_runs)):
        raise ManifestError("run accessions overlap between ENA snapshots")

    runtime: list[dict[str, str]] = []
    metadata_fields = (
        "sample_accession",
        "library_layout",
        "instrument_platform",
        "library_strategy",
        "base_count",
        "read_count",
    )
    for frozen in selected:
        run = frozen["run_accession"]
        ena = indexes[frozen["arm"]].get(run)
        if ena is None:
            raise ManifestError(f"{frozen['arm']}: frozen run {run} absent from snapshot")
        for field in metadata_fields:
            if ena[field] != frozen[field]:
                raise ManifestError(
                    f"{run}: frozen {field}={frozen[field]!r} differs from "
                    f"ENA snapshot {ena[field]!r}"
                )
        urls = _split_pair(ena["fastq_ftp"], run=run, field="fastq_ftp")
        md5s = _split_pair(ena["fastq_md5"], run=run, field="fastq_md5")
        sizes = _split_pair(ena["fastq_bytes"], run=run, field="fastq_bytes")
        if not all(MD5_PATTERN.fullmatch(value) for value in md5s):
            raise ManifestError(f"{run}: invalid FASTQ MD5")
        runtime.append(
            {
                "sample_id": run,
                "cohort": frozen["arm"],
                "accession": run,
                "layout": "PAIRED",
                "fastq1_url": _url(urls[0], run=run),
                "fastq1_md5": md5s[0].lower(),
                "fastq1_bytes": _positive_size(sizes[0], run=run),
                "fastq2_url": _url(urls[1], run=run),
                "fastq2_md5": md5s[1].lower(),
                "fastq2_bytes": _positive_size(sizes[1], run=run),
                "status": "primary",
                "rank": frozen["rank"],
            }
        )
    expected_order = [row["run_accession"] for row in selected]
    if [row["accession"] for row in runtime] != expected_order:
        raise AssertionError("runtime manifest order diverged from frozen selection")
    return runtime


def write_runtime_manifest(path: Path, rows: Sequence[Mapping[str, str]]) -> str:
    payload = canonical_tsv(rows, RUNTIME_MANIFEST_FIELDS)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return hashlib.sha256(payload).hexdigest()


def write_runtime_config(
    *,
    template: Path,
    output: Path,
    manifest: Path,
    project_root: Path,
    scope: str,
) -> str:
    config = yaml.safe_load(template.read_text(encoding="utf-8"))
    try:
        manifest_value = manifest.resolve().relative_to(project_root.resolve()).as_posix()
    except ValueError:
        manifest_value = str(manifest.resolve())
    config["paths"]["sample_manifest"] = manifest_value
    config["experiment"]["id"] = f"{config['experiment']['id']}-{scope}"
    payload = yaml.safe_dump(config, sort_keys=False).encode("utf-8")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(payload)
    return hashlib.sha256(payload).hexdigest()


def prepare_runtime(
    *,
    frozen_manifest: Path,
    snapshot_paths: Mapping[str, Path],
    config_template: Path,
    runtime_manifest: Path,
    runtime_config: Path,
    evidence_output: Path,
    project_root: Path,
    scope: str,
    replacement_ledger: Path | None = None,
) -> dict[str, object]:
    frozen_rows = read_tsv(frozen_manifest)
    snapshots = {arm: read_tsv(path) for arm, path in snapshot_paths.items()}
    replacements: Sequence[Mapping[str, object]] = ()
    if replacement_ledger is not None:
        ledger = yaml.safe_load(replacement_ledger.read_text(encoding="utf-8"))
        if not isinstance(ledger, dict) or ledger.get("schema_version") != 1:
            raise ManifestError("replacement ledger must have schema_version 1")
        candidate = ledger.get("replacements", [])
        if not isinstance(candidate, list) or not all(isinstance(row, dict) for row in candidate):
            raise ManifestError("replacement ledger replacements must be a list")
        replacements = candidate
    runtime_rows = build_runtime_rows(
        frozen_rows, snapshots, scope=scope, replacements=replacements
    )
    manifest_sha = write_runtime_manifest(runtime_manifest, runtime_rows)
    config_sha = write_runtime_config(
        template=config_template,
        output=runtime_config,
        manifest=runtime_manifest,
        project_root=project_root,
        scope=scope,
    )
    per_arm = {}
    for arm in sorted({row["cohort"] for row in runtime_rows}):
        arm_rows = [row for row in runtime_rows if row["cohort"] == arm]
        per_arm[arm] = {
            "runs": len(arm_rows),
            "expected_fastq_bytes": sum(
                int(row["fastq1_bytes"]) + int(row["fastq2_bytes"])
                for row in arm_rows
            ),
        }
    evidence: dict[str, object] = {
        "schema_version": 1,
        "scope": scope,
        "valid": True,
        "ordered_accessions": [row["accession"] for row in runtime_rows],
        "arms": per_arm,
        "inputs": {
            "frozen_manifest_sha256": sha256_bytes(
                canonical_tsv(frozen_rows, MANIFEST_FIELDS)
            ),
            "ena_snapshot_sha256": {
                arm: _sha256(path) for arm, path in sorted(snapshot_paths.items())
            },
            "replacement_ledger_sha256": (
                _sha256(replacement_ledger) if replacement_ledger is not None else None
            ),
        },
        "outputs": {
            "runtime_manifest_sha256": manifest_sha,
            "runtime_config_sha256": config_sha,
        },
        "privacy": {
            "contains_fastq_urls": False,
            "contains_fastq_checksums": False,
            "contains_absolute_paths": False,
            "contains_sequence_data": False,
        },
    }
    write_json_atomic(evidence, evidence_output)
    return evidence
