"""Acquire, verify, union, and index the competitive human reference panel."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import subprocess
import time
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import BinaryIO, Mapping, Protocol, Sequence

from .data_manifest import ManifestError, canonical_tsv, read_tsv, sha256_bytes
from .provenance import write_json_atomic
from .reference_catalog import (
    DONOR_FIELDS,
    PANEL_FIELDS,
    read_csv,
    validate_panel_against_assemblies,
    verify_frozen_selection,
)

METADATA_FIELDS = (
    "source_id",
    "authority",
    "url",
    "expected_bytes",
    "sha256",
    "retrieved_date_utc",
    "fields_used",
)
HEX_MD5 = re.compile(r"^[0-9a-f]{32}$")
HEX_SHA256 = re.compile(r"^[0-9a-f]{64}$")


class AcquisitionError(RuntimeError):
    """An operational reference acquisition or index error."""


class HttpClient(Protocol):
    def text(self, url: str) -> str: ...

    def content_length(self, url: str) -> int: ...

    def download(self, url: str, destination: Path, expected_bytes: int) -> None: ...


class IndexBuilder(Protocol):
    def version(self) -> str: ...

    def build(
        self,
        union_fasta: Path,
        output_index: Path,
        *,
        threads: int,
        index_batch: str,
    ) -> dict[str, object]: ...


class UrlLibClient:
    """HTTPS client with bounded retries and Range-based partial-file resume."""

    def __init__(self, *, timeout_seconds: int = 120, retries: int = 3) -> None:
        self.timeout_seconds = timeout_seconds
        self.retries = retries

    @staticmethod
    def _request(url: str, *, method: str = "GET", headers: Mapping[str, str] | None = None):
        request_headers = {"User-Agent": "hostbias-reference-acquisition/1"}
        request_headers.update(headers or {})
        return urllib.request.Request(url, method=method, headers=request_headers)

    def text(self, url: str) -> str:
        with urllib.request.urlopen(
            self._request(url), timeout=self.timeout_seconds
        ) as response:
            return response.read().decode("utf-8")

    def content_length(self, url: str) -> int:
        with urllib.request.urlopen(
            self._request(url, method="HEAD"), timeout=self.timeout_seconds
        ) as response:
            value = response.headers.get("Content-Length")
        if value is None:
            raise AcquisitionError(f"remote source has no Content-Length: {url}")
        return int(value)

    def download(self, url: str, destination: Path, expected_bytes: int) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)
        partial = destination.with_suffix(destination.suffix + ".partial")
        if partial.exists() and partial.stat().st_size > expected_bytes:
            partial.unlink()
        last_error: BaseException | None = None
        for attempt in range(1, self.retries + 1):
            offset = partial.stat().st_size if partial.exists() else 0
            headers = {"Range": f"bytes={offset}-"} if offset else {}
            try:
                with urllib.request.urlopen(
                    self._request(url, headers=headers),
                    timeout=self.timeout_seconds,
                ) as response:
                    append = offset > 0 and response.status == 206
                    mode = "ab" if append else "wb"
                    with partial.open(mode) as output:
                        while chunk := response.read(8 * 1024 * 1024):
                            output.write(chunk)
                        output.flush()
                        os.fsync(output.fileno())
                if partial.stat().st_size == expected_bytes:
                    os.replace(partial, destination)
                    destination.chmod(0o600)
                    return
                last_error = AcquisitionError(
                    f"downloaded byte count mismatch for {destination.name}"
                )
            except (OSError, urllib.error.URLError) as error:
                last_error = error
            if attempt < self.retries:
                time.sleep(attempt)
        raise AcquisitionError(
            f"download failed after {self.retries} attempts: {destination.name}"
        ) from last_error


def _tool_version(executable: str) -> str:
    try:
        completed = subprocess.run(
            [executable, "--version"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise AcquisitionError(f"cannot execute {executable}") from error
    output = completed.stdout.decode("utf-8", errors="replace")
    line = next((line.strip() for line in output.splitlines() if line.strip()), "")
    if completed.returncode or not line:
        raise AcquisitionError(f"{executable} --version failed")
    return line[:300]


class Minimap2IndexBuilder:
    def __init__(self, executable: str = "minimap2") -> None:
        self.executable = executable

    def version(self) -> str:
        return _tool_version(self.executable)

    def build(
        self,
        union_fasta: Path,
        output_index: Path,
        *,
        threads: int,
        index_batch: str,
    ) -> dict[str, object]:
        output_index.parent.mkdir(parents=True, exist_ok=True)
        partial = output_index.with_suffix(output_index.suffix + ".partial")
        partial.unlink(missing_ok=True)
        stderr_path = output_index.with_suffix(".minimap2.stderr")
        command = [
            self.executable,
            "-x",
            "asm5",
            "-I",
            index_batch,
            "-t",
            str(threads),
            "-d",
            str(partial),
            str(union_fasta),
        ]
        started = time.monotonic()
        with stderr_path.open("wb") as stderr:
            try:
                completed = subprocess.run(
                    command,
                    check=False,
                    stdout=subprocess.DEVNULL,
                    stderr=stderr,
                )
            except OSError as error:
                raise AcquisitionError("minimap2 could not start") from error
        if completed.returncode or not partial.is_file():
            raise AcquisitionError(
                f"minimap2 index build failed with code {completed.returncode}"
            )
        os.replace(partial, output_index)
        output_index.chmod(0o600)
        stderr_path.unlink(missing_ok=True)
        return {
            "preset": "asm5",
            "index_batch": index_batch,
            "threads": threads,
            "seconds": round(time.monotonic() - started, 3),
        }


def multi_hash_file(path: Path) -> dict[str, str]:
    md5 = hashlib.md5(usedforsecurity=False)
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            md5.update(chunk)
            sha256.update(chunk)
    return {"md5": md5.hexdigest(), "sha256": sha256.hexdigest()}


def _expected_md5_from_sidecar(text: str, filename: str) -> str:
    matches = []
    for line in text.splitlines():
        fields = line.strip().split()
        if len(fields) >= 2 and Path(fields[-1].lstrip("*")).name == filename:
            digest = fields[0].lower()
            if HEX_MD5.fullmatch(digest):
                matches.append(digest)
    if len(matches) != 1:
        raise AcquisitionError(
            f"expected one authoritative MD5 for {filename}; found {len(matches)}"
        )
    return matches[0]


def _verify_existing_reference(
    path: Path, row: Mapping[str, str]
) -> dict[str, object] | None:
    if not path.is_file() or path.stat().st_size != int(row["expected_bytes"]):
        return None
    hashes = multi_hash_file(path)
    if hashes["md5"] != row["expected_md5"]:
        return None
    expected_sha256 = row["expected_sha256"]
    if expected_sha256 and hashes["sha256"] != expected_sha256:
        return None
    return {
        "reference_id": row["reference_id"],
        "kind": row["kind"],
        "donor_id": row["donor_id"],
        "haplotype": row["haplotype"],
        "superpopulation": row["superpopulation"],
        "panel_role": row["panel_role"],
        "genbank_accession": row["genbank_accession"],
        "bytes": path.stat().st_size,
        "md5": hashes["md5"],
        "sha256": hashes["sha256"],
        "include_in_union": row["include_in_union"] == "true",
        "filename": path.name,
    }


def acquire_reference(
    row: Mapping[str, str],
    destination_root: Path,
    client: HttpClient,
) -> dict[str, object]:
    destination = destination_root / row["local_filename"]
    expected_bytes = int(row["expected_bytes"])
    existing = _verify_existing_reference(destination, row)
    if existing is not None:
        existing["reused"] = True
        return existing
    if destination.exists():
        destination.unlink()
    remote_bytes = client.content_length(row["source_url"])
    if remote_bytes != expected_bytes:
        raise AcquisitionError(
            f"{row['reference_id']}: remote bytes changed "
            f"({remote_bytes} != {expected_bytes})"
        )
    remote_md5 = _expected_md5_from_sidecar(
        client.text(row["source_md5_url"]), row["local_filename"]
    )
    if remote_md5 != row["expected_md5"]:
        raise AcquisitionError(f"{row['reference_id']}: authoritative MD5 changed")
    client.download(row["source_url"], destination, expected_bytes)
    verified = _verify_existing_reference(destination, row)
    if verified is None:
        destination.unlink(missing_ok=True)
        raise AcquisitionError(f"{row['reference_id']}: local checksum verification failed")
    verified["reused"] = False
    return verified


def acquire_metadata_sources(
    rows: Sequence[Mapping[str, str]],
    metadata_root: Path,
    client: HttpClient,
) -> tuple[list[dict[str, object]], dict[str, Path]]:
    metadata_root.mkdir(parents=True, exist_ok=True)
    evidence = []
    paths: dict[str, Path] = {}
    for row in rows:
        expected_bytes = int(row["expected_bytes"])
        destination = metadata_root / f"{row['source_id']}.metadata"
        valid = (
            destination.is_file()
            and destination.stat().st_size == expected_bytes
            and multi_hash_file(destination)["sha256"] == row["sha256"]
        )
        if not valid:
            destination.unlink(missing_ok=True)
            client.download(row["url"], destination, expected_bytes)
        observed = multi_hash_file(destination)
        if (
            destination.stat().st_size != expected_bytes
            or observed["sha256"] != row["sha256"]
        ):
            raise AcquisitionError(f"{row['source_id']}: metadata verification failed")
        paths[row["source_id"]] = destination
        evidence.append(
            {
                "source_id": row["source_id"],
                "authority": row["authority"],
                "url": row["url"],
                "bytes": expected_bytes,
                "sha256": observed["sha256"],
            }
        )
    return evidence, paths


class _HashingWriter:
    def __init__(self, handle: BinaryIO) -> None:
        self.handle = handle
        self.digest = hashlib.sha256()
        self.bytes = 0

    def write(self, payload: bytes) -> None:
        self.handle.write(payload)
        self.digest.update(payload)
        self.bytes += len(payload)


def build_union_fasta(
    panel: Sequence[Mapping[str, str]],
    downloads_root: Path,
    output_fasta: Path,
) -> dict[str, object]:
    """Concatenate primary references while namespacing every FASTA header."""
    included = [row for row in panel if row["include_in_union"] == "true"]
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    partial = output_fasta.with_suffix(output_fasta.suffix + ".partial")
    partial.unlink(missing_ok=True)
    sequence_count = 0
    base_count = 0
    headers: set[bytes] = set()
    with partial.open("wb") as raw_output:
        writer = _HashingWriter(raw_output)
        for row in included:
            source = downloads_root / row["local_filename"]
            with gzip.open(source, "rb") as fasta:
                saw_header = False
                for line in fasta:
                    if line.startswith(b">"):
                        saw_header = True
                        original = line[1:].rstrip(b"\r\n")
                        if not original:
                            raise AcquisitionError(
                                f"{row['reference_id']}: empty FASTA header"
                            )
                        header = row["reference_id"].encode("ascii") + b"|" + original
                        token = header.split(maxsplit=1)[0]
                        if token in headers:
                            raise AcquisitionError(
                                f"{row['reference_id']}: duplicate FASTA header"
                            )
                        headers.add(token)
                        writer.write(b">" + header + b"\n")
                        sequence_count += 1
                    else:
                        if not saw_header:
                            raise AcquisitionError(
                                f"{row['reference_id']}: sequence precedes FASTA header"
                            )
                        sequence = line.strip()
                        base_count += len(sequence)
                        writer.write(sequence + b"\n")
            if not saw_header:
                raise AcquisitionError(f"{row['reference_id']}: empty FASTA")
        raw_output.flush()
        os.fsync(raw_output.fileno())
    os.replace(partial, output_fasta)
    output_fasta.chmod(0o600)
    return {
        "reference_id": "hprc-balanced-chm13-v1",
        "filename": output_fasta.name,
        "included_references": [row["reference_id"] for row in included],
        "included_reference_count": len(included),
        "sequence_count": sequence_count,
        "base_count": base_count,
        "bytes": writer.bytes,
        "sha256": writer.digest.hexdigest(),
        "headers_namespaced": True,
    }


def _load_json(path: Path) -> dict[str, object] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def build_reference_panel(
    *,
    metadata_sources_path: Path,
    donors_path: Path,
    panel_path: Path,
    reference_root: Path,
    checkpoint_path: Path,
    threads: int,
    index_batch: str = "64G",
    client: HttpClient | None = None,
    index_builder: IndexBuilder | None = None,
) -> dict[str, object]:
    """Execute or resume the complete reference acquisition and index build."""
    client = client or UrlLibClient()
    index_builder = index_builder or Minimap2IndexBuilder()
    metadata_sources = read_tsv(metadata_sources_path)
    donors = read_tsv(donors_path)
    panel = read_tsv(panel_path)
    if set(METADATA_FIELDS).difference(metadata_sources[0]):
        raise ManifestError("reference metadata source ledger has missing fields")

    metadata_evidence, metadata_paths = acquire_metadata_sources(
        metadata_sources, reference_root / "metadata", client
    )
    selection = verify_frozen_selection(
        assembly_metadata=metadata_paths["hprc_release2_assemblies"],
        sample_metadata=metadata_paths["hprc_release2_samples"],
        igsr_panel=metadata_paths["igsr_1000g_20130502_panel"],
        frozen_donors=donors_path,
        panel_manifest=panel_path,
    )
    assemblies = read_csv(metadata_paths["hprc_release2_assemblies"])
    validate_panel_against_assemblies(panel, assemblies)

    downloads_root = reference_root / "downloads"
    downloads_root.mkdir(parents=True, exist_ok=True)
    download_evidence = [
        acquire_reference(row, downloads_root, client) for row in panel
    ]
    input_fingerprint = sha256_bytes(
        canonical_tsv(panel, PANEL_FIELDS)
        + canonical_tsv(donors, DONOR_FIELDS)
        + "\n".join(
            f"{row['reference_id']}\t{row['sha256']}" for row in download_evidence
        ).encode("ascii")
    )

    union_fasta = reference_root / "panel" / "hprc-balanced-chm13-v1.fa"
    union_checkpoint_path = reference_root / "panel" / "union.checkpoint.json"
    union_checkpoint = _load_json(union_checkpoint_path)
    if (
        union_checkpoint is None
        or union_checkpoint.get("input_fingerprint") != input_fingerprint
        or not union_fasta.is_file()
        or union_fasta.stat().st_size != union_checkpoint.get("bytes")
        or multi_hash_file(union_fasta)["sha256"] != union_checkpoint.get("sha256")
    ):
        union_evidence = build_union_fasta(panel, downloads_root, union_fasta)
        union_checkpoint = {**union_evidence, "input_fingerprint": input_fingerprint}
        write_json_atomic(union_checkpoint, union_checkpoint_path)
    union_evidence = {
        key: value
        for key, value in union_checkpoint.items()
        if key != "input_fingerprint"
    }

    minimap2_version = index_builder.version()
    index_path = reference_root / "indexes" / "hprc-balanced-chm13-v1.mmi"
    index_checkpoint_path = reference_root / "indexes" / "index.checkpoint.json"
    index_checkpoint = _load_json(index_checkpoint_path)
    index_reusable = (
        index_checkpoint is not None
        and index_checkpoint.get("union_sha256") == union_evidence["sha256"]
        and index_checkpoint.get("minimap2_version") == minimap2_version
        and index_checkpoint.get("index_batch") == index_batch
        and index_path.is_file()
        and index_path.stat().st_size == index_checkpoint.get("bytes")
        and multi_hash_file(index_path)["sha256"] == index_checkpoint.get("sha256")
    )
    if not index_reusable:
        build_details = index_builder.build(
            union_fasta,
            index_path,
            threads=threads,
            index_batch=index_batch,
        )
        index_hash = multi_hash_file(index_path)["sha256"]
        index_checkpoint = {
            "reference_id": "hprc-balanced-chm13-v1",
            "filename": index_path.name,
            "bytes": index_path.stat().st_size,
            "sha256": index_hash,
            "union_sha256": union_evidence["sha256"],
            "minimap2_version": minimap2_version,
            **build_details,
        }
        write_json_atomic(index_checkpoint, index_checkpoint_path)

    report: dict[str, object] = {
        "schema_version": 1,
        "checkpoint": "competitive_human_reference_panel",
        "captured_at_utc": datetime.now(UTC).isoformat(),
        "status": "PASS",
        "reference_id": "hprc-balanced-chm13-v1",
        "provenance": {
            "metadata_sources_manifest_sha256": sha256_bytes(
                canonical_tsv(metadata_sources, METADATA_FIELDS)
            ),
            "donor_catalog_sha256": sha256_bytes(canonical_tsv(donors, DONOR_FIELDS)),
            "panel_manifest_sha256": sha256_bytes(canonical_tsv(panel, PANEL_FIELDS)),
            "selection": selection,
            "metadata_sources": metadata_evidence,
        },
        "downloads": download_evidence,
        "download_totals": {
            "references": len(download_evidence),
            "bytes": sum(int(row["bytes"]) for row in download_evidence),
            "reused": sum(bool(row["reused"]) for row in download_evidence),
        },
        "union_fasta": union_evidence,
        "minimap2_index": index_checkpoint,
        "privacy": {
            "reference_sequences_committed": False,
            "absolute_paths_recorded": False,
            "aggregate_evidence_only": True,
        },
    }
    write_json_atomic(report, checkpoint_path)
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-sources", required=True, type=Path)
    parser.add_argument("--donors", required=True, type=Path)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--index-batch", default="64G")
    parser.add_argument("--minimap2", default="minimap2")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    report = build_reference_panel(
        metadata_sources_path=args.metadata_sources,
        donors_path=args.donors,
        panel_path=args.panel,
        reference_root=args.reference_root,
        checkpoint_path=args.checkpoint,
        threads=args.threads,
        index_batch=args.index_batch,
        index_builder=Minimap2IndexBuilder(args.minimap2),
    )
    print(
        json.dumps(
            {
                "status": report["status"],
                "reference_id": report["reference_id"],
                "checkpoint": str(args.checkpoint),
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
