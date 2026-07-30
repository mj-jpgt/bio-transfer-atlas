"""Convert real minimap2 PAF results into the competitive-label input contract."""

from __future__ import annotations

import csv
import json
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from hostbias.assembly_qc import iter_fasta
from hostbias.provenance import sha256_file, write_json_atomic
from hostbias.schemas import AlignmentRow, SchemaError


@dataclass(frozen=True)
class AlignmentSource:
    domain: str
    reference_id: str
    reference_sha256: str
    paf_path: Path


@dataclass(frozen=True)
class AlignmentRunSpec:
    sample_id: str
    filter_mode: str
    assembly_path: Path
    assembly_sha256: str
    mapper_version: str
    mapper_preset: str
    sources: tuple[AlignmentSource, ...]


def _required_text(mapping: dict[str, Any], key: str, context: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value:
        raise SchemaError(f"{context}.{key} must be a non-empty string")
    return value


def load_alignment_run_spec(path: Path) -> AlignmentRunSpec:
    """Load the private mapper run specification and resolve its relative paths."""

    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SchemaError(f"cannot read alignment run spec {path}: {exc}") from exc
    if not isinstance(raw, dict) or raw.get("schema_version") != "1.0":
        raise SchemaError("alignment run spec schema_version must be '1.0'")
    expected = {
        "schema_version",
        "sample_id",
        "filter_mode",
        "assembly",
        "mapper",
        "references",
    }
    if set(raw) != expected:
        raise SchemaError(
            f"alignment run spec keys must be exactly {sorted(expected)}"
        )
    filter_mode = _required_text(raw, "filter_mode", "spec")
    if filter_mode not in {"source", "strict"}:
        raise SchemaError("spec.filter_mode must be 'source' or 'strict'")
    assembly = raw["assembly"]
    mapper = raw["mapper"]
    references = raw["references"]
    if not isinstance(assembly, dict) or not isinstance(mapper, dict):
        raise SchemaError("spec assembly and mapper must be objects")
    if not isinstance(references, list) or not references:
        raise SchemaError("spec.references must be a non-empty list")
    assembly_path = (path.parent / _required_text(assembly, "path", "assembly")).resolve()
    sources: list[AlignmentSource] = []
    for index, reference in enumerate(references):
        context = f"references[{index}]"
        if not isinstance(reference, dict):
            raise SchemaError(f"{context} must be an object")
        domain = _required_text(reference, "domain", context)
        if domain not in {"human", "gtdb"}:
            raise SchemaError(f"{context}.domain must be 'human' or 'gtdb'")
        sources.append(
            AlignmentSource(
                domain=domain,
                reference_id=_required_text(reference, "reference_id", context),
                reference_sha256=_required_text(
                    reference, "reference_sha256", context
                ),
                paf_path=(
                    path.parent / _required_text(reference, "paf_path", context)
                ).resolve(),
            )
        )
    domains = {source.domain for source in sources}
    if domains != {"human", "gtdb"}:
        raise SchemaError("spec.references must include human and gtdb sources")
    spec = AlignmentRunSpec(
        sample_id=_required_text(raw, "sample_id", "spec"),
        filter_mode=filter_mode,
        assembly_path=assembly_path,
        assembly_sha256=_required_text(assembly, "sha256", "assembly"),
        mapper_version=_required_text(mapper, "version", "mapper"),
        mapper_preset=_required_text(mapper, "preset", "mapper"),
        sources=tuple(sources),
    )
    actual_assembly_sha256 = sha256_file(spec.assembly_path)
    if actual_assembly_sha256 != spec.assembly_sha256:
        raise SchemaError(
            "assembly checksum mismatch: "
            f"expected {spec.assembly_sha256}, found {actual_assembly_sha256}"
        )
    for source in spec.sources:
        if len(source.reference_sha256) != 64 or any(
            character not in "0123456789abcdef"
            for character in source.reference_sha256
        ):
            raise SchemaError(
                f"invalid reference sha256 for {source.reference_id!r}"
            )
    return spec


def _parse_score(tags: list[str], path: Path, line_number: int) -> float:
    for tag in tags:
        if tag.startswith("AS:i:"):
            score = float(tag.removeprefix("AS:i:"))
            if math.isfinite(score) and score >= 0:
                return score
            break
    raise SchemaError(f"{path}:{line_number}: missing/invalid non-negative AS:i tag")


def _parse_paf(
    source: AlignmentSource,
    sample_id: str,
    assembly_lengths: dict[str, int],
) -> list[AlignmentRow]:
    rows: list[AlignmentRow] = []
    with source.paf_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise SchemaError(
                    f"{source.paf_path}:{line_number}: PAF requires >= 12 columns"
                )
            query_id = fields[0]
            if query_id not in assembly_lengths:
                raise SchemaError(
                    f"{source.paf_path}:{line_number}: unknown query {query_id!r}"
                )
            try:
                query_length = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                matching_bases = int(fields[9])
                alignment_block_length = int(fields[10])
                mapq = float(fields[11])
            except ValueError as exc:
                raise SchemaError(
                    f"{source.paf_path}:{line_number}: non-numeric PAF field"
                ) from exc
            if query_length != assembly_lengths[query_id]:
                raise SchemaError(
                    f"{source.paf_path}:{line_number}: query length mismatch for "
                    f"{query_id!r}"
                )
            if not (
                0 <= query_start < query_end <= query_length
                and 0 <= matching_bases <= alignment_block_length
                and alignment_block_length > 0
            ):
                raise SchemaError(
                    f"{source.paf_path}:{line_number}: invalid PAF coordinates"
                )
            aligned_query_bp = query_end - query_start
            rows.append(
                AlignmentRow(
                    sample_id=sample_id,
                    contig_id=query_id,
                    contig_length=query_length,
                    target_domain=source.domain,
                    aligned_bp=aligned_query_bp,
                    identity=matching_bases / alignment_block_length,
                    query_coverage=aligned_query_bp / query_length,
                    mapq=mapq,
                    alignment_score=_parse_score(
                        fields[12:], source.paf_path, line_number
                    ),
                )
            )
    return rows


def build_alignment_table(
    spec: AlignmentRunSpec,
) -> tuple[list[AlignmentRow], dict[str, Any]]:
    """Build sensitive contig rows plus a separately publishable run manifest."""

    assembly_lengths = {
        identifier: len(sequence) for identifier, sequence in iter_fasta(spec.assembly_path)
    }
    rows: list[AlignmentRow] = []
    hit_queries: dict[str, set[str]] = {"human": set(), "gtdb": set()}
    source_summaries: list[dict[str, Any]] = []
    for source in spec.sources:
        source_rows = _parse_paf(source, spec.sample_id, assembly_lengths)
        rows.extend(source_rows)
        hit_queries[source.domain].update(row.contig_id for row in source_rows)
        source_summaries.append(
            {
                "domain": source.domain,
                "reference_id": source.reference_id,
                "reference_sha256": source.reference_sha256,
                "paf_sha256": sha256_file(source.paf_path),
                "alignment_row_count": len(source_rows),
                "query_with_hit_count": len(
                    {row.contig_id for row in source_rows}
                ),
            }
        )
    any_hits = hit_queries["human"] | hit_queries["gtdb"]
    for query_id, query_length in assembly_lengths.items():
        if query_id not in any_hits:
            rows.append(
                AlignmentRow(
                    spec.sample_id,
                    query_id,
                    query_length,
                    "none",
                    0,
                    0,
                    0,
                    0,
                    0,
                )
            )
    rows.sort(
        key=lambda row: (
            row.sample_id,
            row.contig_id,
            row.target_domain,
            -row.alignment_score,
        )
    )
    manifest = {
        "schema_version": "1.0",
        "sample_id": spec.sample_id,
        "filter_mode": spec.filter_mode,
        "assembly_sha256": spec.assembly_sha256,
        "assembly_contig_count": len(assembly_lengths),
        "mapper": {
            "name": "minimap2",
            "version": spec.mapper_version,
            "preset": spec.mapper_preset,
        },
        "references": sorted(
            source_summaries,
            key=lambda row: (row["domain"], row["reference_id"]),
        ),
        "aggregate": {
            "alignment_row_count": len(rows),
            "human_hit_contig_count": len(hit_queries["human"]),
            "gtdb_hit_contig_count": len(hit_queries["gtdb"]),
            "no_hit_contig_count": len(assembly_lengths) - len(any_hits),
        },
        "privacy": {
            "contains_sequences": False,
            "contains_contig_identifiers": False,
            "contains_filesystem_paths": False,
        },
    }
    return rows, manifest


def write_alignment_tsv_atomic(rows: list[AlignmentRow], output: Path) -> None:
    """Write the sensitive analysis handoff atomically."""

    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp"
    )
    header = list(AlignmentRow.PARSERS)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(
                    {
                        "sample_id": row.sample_id,
                        "contig_id": row.contig_id,
                        "contig_length": row.contig_length,
                        "target_domain": row.target_domain,
                        "aligned_bp": row.aligned_bp,
                        "identity": f"{row.identity:.12g}",
                        "query_coverage": f"{row.query_coverage:.12g}",
                        "mapq": f"{row.mapq:.12g}",
                        "alignment_score": f"{row.alignment_score:.12g}",
                    }
                )
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, output)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def run_alignment_bridge(
    spec_path: Path, output_tsv: Path, output_manifest: Path
) -> None:
    spec = load_alignment_run_spec(spec_path)
    rows, manifest = build_alignment_table(spec)
    write_alignment_tsv_atomic(rows, output_tsv)
    write_json_atomic(manifest, output_manifest)
