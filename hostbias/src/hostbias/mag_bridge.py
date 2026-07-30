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


def _normalize_bin_id(value: str) -> str:
    value = Path(value.strip()).name
    for suffix in (".fasta", ".fna", ".fa"):
        if value.lower().endswith(suffix):
            return value[: -len(suffix)]
    return value


def _read_dict_tsv(
    path: Path, *, allow_empty: bool = False
) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SchemaError(f"{path}: missing TSV header")
        rows = [dict(row) for row in reader]
    if not rows and not allow_empty:
        raise SchemaError(f"{path}: table has no data rows")
    return rows


def _parse_dastool_map(path: Path) -> list[tuple[str, str]]:
    assignments: list[tuple[str, str]] = []
    seen_contigs: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if row == ["contig_id", "bin_id"] or row == ["scaffold", "bin"]:
                continue
            if len(row) != 2 or not all(row):
                raise SchemaError(
                    f"{path}:{line_number}: expected two non-empty columns"
                )
            contig_id, bin_id = row[0], _normalize_bin_id(row[1])
            if contig_id in seen_contigs:
                raise SchemaError(
                    f"{path}:{line_number}: duplicate selected contig {contig_id!r}"
                )
            seen_contigs.add(contig_id)
            assignments.append((contig_id, bin_id))
    return sorted(assignments)


def _parse_checkm2(path: Path) -> dict[str, tuple[float, float]]:
    results: dict[str, tuple[float, float]] = {}
    for line_number, row in enumerate(_read_dict_tsv(path), start=2):
        try:
            bin_id = _normalize_bin_id(row["Name"])
            completeness = float(row["Completeness"])
            contamination = float(row["Contamination"])
        except (KeyError, ValueError) as exc:
            raise SchemaError(
                f"{path}:{line_number}: invalid CheckM2 row"
            ) from exc
        if (
            not bin_id
            or bin_id in results
            or not 0 <= completeness <= 100
            or contamination < 0
        ):
            raise SchemaError(
                f"{path}:{line_number}: invalid/duplicate CheckM2 result"
            )
        results[bin_id] = (completeness / 100, contamination / 100)
    return results


def _parse_gunc(path: Path) -> dict[str, bool]:
    results: dict[str, bool] = {}
    for line_number, row in enumerate(_read_dict_tsv(path), start=2):
        try:
            bin_id = _normalize_bin_id(row["genome"])
            raw_pass = row["pass.GUNC"].strip().lower()
        except KeyError as exc:
            raise SchemaError(f"{path}:{line_number}: invalid GUNC row") from exc
        if raw_pass in {"true", "1"}:
            passed = True
        elif raw_pass in {"false", "0", "nan", "na", ""}:
            passed = False
        else:
            raise SchemaError(
                f"{path}:{line_number}: invalid pass.GUNC value {raw_pass!r}"
            )
        if not bin_id or bin_id in results:
            raise SchemaError(f"{path}:{line_number}: duplicate/empty GUNC genome")
        results[bin_id] = passed
    return results


def _taxonomy_ranks(classification: str) -> tuple[str, str | None, str | None]:
    ranks = {
        token[:3]: token
        for token in classification.strip().split(";")
        if len(token) >= 3 and token[1:3] == "__"
    }
    domain_token = ranks.get("d__")
    domain = (
        domain_token.removeprefix("d__")
        if domain_token in {"d__Bacteria", "d__Archaea"}
        else "Unclassified"
    )
    genus = ranks.get("g__")
    species = ranks.get("s__")
    if genus == "g__":
        genus = None
    if species == "s__":
        species = None
    return domain, genus, species


def _parse_gtdb(paths: Iterable[Path]) -> dict[str, tuple[str, str | None, str | None]]:
    results: dict[str, tuple[str, str | None, str | None]] = {}
    for path in paths:
        for line_number, row in enumerate(
            _read_dict_tsv(path, allow_empty=True), start=2
        ):
            try:
                bin_id = _normalize_bin_id(row["user_genome"])
                taxonomy = _taxonomy_ranks(row["classification"])
            except KeyError as exc:
                raise SchemaError(
                    f"{path}:{line_number}: invalid GTDB-Tk summary row"
                ) from exc
            if not bin_id or bin_id in results:
                raise SchemaError(
                    f"{path}:{line_number}: duplicate/empty GTDB-Tk genome"
                )
            results[bin_id] = taxonomy
    return results


def build_mag_contracts(
    sample_id: str,
    dastool_map: Path,
    checkm2_report: Path,
    gunc_report: Path,
    gtdb_summaries: tuple[Path, ...],
    contig_bins_output: Path,
    bin_qc_output: Path,
) -> int:
    """Join DAS Tool, CheckM2, GUNC, and GTDB-Tk with exact bin coverage."""

    assignments = _parse_dastool_map(dastool_map)
    selected_bins = {bin_id for _, bin_id in assignments}
    if not selected_bins:
        _atomic_tsv(
            contig_bins_output,
            ("sample_id", "contig_id", "bin_id"),
            (),
        )
        _atomic_tsv(
            bin_qc_output,
            (
                "sample_id",
                "bin_id",
                "das_tool_selected",
                "checkm2_completeness",
                "checkm2_contamination",
                "gunc_pass",
                "gtdb_domain",
                "gtdb_genus",
                "gtdb_species",
            ),
            (),
        )
        return 0
    checkm2 = _parse_checkm2(checkm2_report)
    gunc = _parse_gunc(gunc_report)
    gtdb = _parse_gtdb(gtdb_summaries)
    for label, observed in (
        ("CheckM2", set(checkm2)),
        ("GUNC", set(gunc)),
        ("GTDB-Tk", set(gtdb)),
    ):
        if observed != selected_bins:
            raise SchemaError(
                f"{label} bins differ from DAS Tool selection; "
                f"missing={sorted(selected_bins - observed)}, "
                f"extra={sorted(observed - selected_bins)}"
            )

    _atomic_tsv(
        contig_bins_output,
        ("sample_id", "contig_id", "bin_id"),
        ((sample_id, contig_id, bin_id) for contig_id, bin_id in assignments),
    )
    qc_rows = []
    for bin_id in sorted(selected_bins):
        completeness, contamination = checkm2[bin_id]
        domain, genus, species = gtdb[bin_id]
        qc_rows.append(
            (
                sample_id,
                bin_id,
                "true",
                f"{completeness:.12g}",
                f"{contamination:.12g}",
                str(gunc[bin_id]).lower(),
                domain,
                genus or "",
                species or "",
            )
        )
    _atomic_tsv(
        bin_qc_output,
        (
            "sample_id",
            "bin_id",
            "das_tool_selected",
            "checkm2_completeness",
            "checkm2_contamination",
            "gunc_pass",
            "gtdb_domain",
            "gtdb_genus",
            "gtdb_species",
        ),
        qc_rows,
    )
    return len(selected_bins)
