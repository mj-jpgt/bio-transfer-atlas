"""Streaming, identifier-free assembly quality-control aggregates."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from typing import Any

from hostbias.provenance import sha256_file
from hostbias.schemas import SchemaError


IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield unique FASTA identifiers and uppercase sequences."""

    identifier: str | None = None
    chunks: list[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(chunks)
                identifier = stripped[1:].split(maxsplit=1)[0]
                if not identifier:
                    raise SchemaError(f"{path}:{line_number}: empty FASTA identifier")
                if identifier in seen:
                    raise SchemaError(
                        f"{path}:{line_number}: duplicate FASTA identifier {identifier!r}"
                    )
                seen.add(identifier)
                chunks = []
            elif identifier is None:
                raise SchemaError(
                    f"{path}:{line_number}: sequence occurs before first header"
                )
            else:
                sequence = stripped.upper()
                invalid = set(sequence) - IUPAC_DNA
                if invalid:
                    raise SchemaError(
                        f"{path}:{line_number}: invalid DNA symbols {sorted(invalid)}"
                    )
                chunks.append(sequence)
    if identifier is not None:
        yield identifier, "".join(chunks)


def _nx(lengths: list[int], fraction: float) -> int:
    threshold = sum(lengths) * fraction
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= threshold:
            return length
    raise AssertionError("Nx requires non-empty lengths")


def assembly_qc(path: Path, sample_id: str, filter_mode: str) -> dict[str, Any]:
    """Return aggregate assembly metrics without sequence or contig identifiers."""

    if filter_mode not in {"source", "strict"}:
        raise SchemaError("filter_mode must be 'source' or 'strict'")
    lengths: list[int] = []
    gc_bases = 0
    acgt_bases = 0
    n_bases = 0
    ambiguous_bases = 0
    for _, sequence in iter_fasta(path):
        if not sequence:
            raise SchemaError(f"{path}: empty contig is not permitted")
        lengths.append(len(sequence))
        gc_bases += sequence.count("G") + sequence.count("C")
        acgt_bases += sum(sequence.count(base) for base in "ACGT")
        n_bases += sequence.count("N")
        ambiguous_bases += sum(sequence.count(base) for base in "RYSWKMBDHV")
    if not lengths:
        raise SchemaError(f"{path}: assembly contains no contigs")
    total_bp = sum(lengths)
    return {
        "schema_version": "1.0",
        "sample_id": sample_id,
        "filter_mode": filter_mode,
        "assembly_sha256": sha256_file(path),
        "contig_count": len(lengths),
        "total_bp": total_bp,
        "minimum_contig_bp": min(lengths),
        "maximum_contig_bp": max(lengths),
        "n50_bp": _nx(lengths, 0.50),
        "n90_bp": _nx(lengths, 0.90),
        "gc_fraction_acgt": gc_bases / acgt_bases if acgt_bases else 0.0,
        "n_fraction": n_bases / total_bp,
        "other_ambiguous_fraction": ambiguous_bases / total_bp,
        "privacy": {
            "contains_sequences": False,
            "contains_contig_identifiers": False,
            "contains_filesystem_paths": False,
        },
    }
