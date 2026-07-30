"""Strict translators between MAG tools and HostBias private TSV contracts."""

from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path
from typing import Iterable

from hostbias.assembly_qc import iter_fasta
from hostbias.schemas import SchemaError


def _atomic_tsv(
    output: Path, header: tuple[str, ...] | None, rows: Iterable[tuple[object, ...]]
) -> None:
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            if header is not None:
                writer.writerow(header)
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, output)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def depth_to_maxbin_abundance(depth_path: Path, output: Path) -> int:
    """Convert MetaBAT's depth table to MaxBin's two-column abundance format."""

    abundance: list[tuple[str, str]] = []
    with depth_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"contigName", "contigLen", "totalAvgDepth"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise SchemaError(
                f"{depth_path}: depth header must contain {sorted(required)}"
            )
        seen: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            contig = row["contigName"]
            if not contig or contig in seen:
                raise SchemaError(
                    f"{depth_path}:{line_number}: empty or duplicate contig"
                )
            seen.add(contig)
            try:
                length = int(row["contigLen"])
                value = float(row["totalAvgDepth"])
            except ValueError as exc:
                raise SchemaError(
                    f"{depth_path}:{line_number}: invalid length/depth"
                ) from exc
            if length <= 0 or value < 0:
                raise SchemaError(
                    f"{depth_path}:{line_number}: length must be positive and "
                    "depth non-negative"
                )
            abundance.append((contig, f"{value:.12g}"))
    if not abundance:
        raise SchemaError(f"{depth_path}: depth table has no contigs")
    _atomic_tsv(output, None, abundance)
    return len(abundance)


def bins_to_scaffolds2bin(
    bin_dir: Path,
    output: Path,
    bin_prefix: str,
    extensions: tuple[str, ...] = (".fa", ".fna", ".fasta"),
) -> int:
    """Convert FASTA bins to DAS Tool's headerless scaffold-to-bin contract."""

    if not bin_prefix or any(character in bin_prefix for character in "\t\r\n"):
        raise SchemaError("bin_prefix must be non-empty and tab/newline-free")
    files = sorted(
        (
            path
            for path in bin_dir.iterdir()
            if path.is_file() and path.suffix.lower() in extensions
        ),
        key=lambda path: path.name,
    )
    assignments: list[tuple[str, str]] = []
    assigned: set[str] = set()
    for path in files:
        bin_id = f"{bin_prefix}{path.stem}"
        for contig_id, _ in iter_fasta(path):
            if contig_id in assigned:
                raise SchemaError(
                    f"contig {contig_id!r} occurs in multiple bins under {bin_dir}"
                )
            assigned.add(contig_id)
            assignments.append((contig_id, bin_id))
    assignments.sort()
    _atomic_tsv(output, None, assignments)
    return len(assignments)
