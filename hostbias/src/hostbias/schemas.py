"""Strict tabular schemas used at the workflow/analysis boundary."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Callable, ClassVar, Iterable, TypeVar


class SchemaError(ValueError):
    """Raised when an input table cannot be interpreted unambiguously."""


def _text(value: str) -> str:
    value = value.strip()
    if not value:
        raise ValueError("must not be empty")
    return value


def _nonnegative_int(value: str) -> int:
    result = int(value)
    if result < 0:
        raise ValueError("must be non-negative")
    return result


def _finite_nonnegative(value: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result < 0:
        raise ValueError("must be finite and non-negative")
    return result


def _fraction(value: str) -> float:
    result = float(value)
    if not math.isfinite(result) or not 0 <= result <= 1:
        raise ValueError("must be in [0, 1]")
    return result


def _boolean(value: str) -> bool:
    lowered = value.strip().lower()
    if lowered in {"true", "1", "yes"}:
        return True
    if lowered in {"false", "0", "no"}:
        return False
    raise ValueError("must be a boolean")


def _optional_text(value: str) -> str | None:
    value = value.strip()
    return value or None


@dataclass(frozen=True)
class AlignmentRow:
    sample_id: str
    contig_id: str
    contig_length: int
    target_domain: str
    aligned_bp: int
    identity: float
    query_coverage: float
    mapq: float
    alignment_score: float

    PARSERS: ClassVar[dict[str, Callable[[str], object]]]

    def __post_init__(self) -> None:
        if self.target_domain not in {"human", "gtdb"}:
            raise SchemaError("target_domain must be 'human' or 'gtdb'")
        if self.contig_length <= 0:
            raise SchemaError("contig_length must be positive")
        if self.aligned_bp > self.contig_length:
            raise SchemaError("aligned_bp cannot exceed contig_length")


AlignmentRow.PARSERS = {
    "sample_id": _text,
    "contig_id": _text,
    "contig_length": _nonnegative_int,
    "target_domain": _text,
    "aligned_bp": _nonnegative_int,
    "identity": _fraction,
    "query_coverage": _fraction,
    "mapq": _finite_nonnegative,
    "alignment_score": _finite_nonnegative,
}


@dataclass(frozen=True)
class ContigBinRow:
    sample_id: str
    contig_id: str
    bin_id: str

    PARSERS: ClassVar[dict[str, Callable[[str], object]]]


ContigBinRow.PARSERS = {
    "sample_id": _text,
    "contig_id": _text,
    "bin_id": _text,
}


@dataclass(frozen=True)
class BinQcRow:
    sample_id: str
    bin_id: str
    das_tool_selected: bool
    checkm2_completeness: float
    checkm2_contamination: float
    gunc_pass: bool
    gtdb_domain: str
    gtdb_genus: str | None
    gtdb_species: str | None

    PARSERS: ClassVar[dict[str, Callable[[str], object]]]

    def __post_init__(self) -> None:
        if self.gtdb_domain not in {"Bacteria", "Archaea", "Eukaryota", "Unclassified"}:
            raise SchemaError(f"unsupported gtdb_domain: {self.gtdb_domain}")


BinQcRow.PARSERS = {
    "sample_id": _text,
    "bin_id": _text,
    "das_tool_selected": _boolean,
    "checkm2_completeness": _fraction,
    "checkm2_contamination": _fraction,
    "gunc_pass": _boolean,
    "gtdb_domain": _text,
    "gtdb_genus": _optional_text,
    "gtdb_species": _optional_text,
}


@dataclass(frozen=True)
class ControlTruthRow:
    sample_id: str
    contig_id: str
    truth: str
    contig_length: int

    PARSERS: ClassVar[dict[str, Callable[[str], object]]]

    def __post_init__(self) -> None:
        if self.truth not in {"human", "microbial"}:
            raise SchemaError("truth must be 'human' or 'microbial'")
        if self.contig_length <= 0:
            raise SchemaError("contig_length must be positive")


ControlTruthRow.PARSERS = {
    "sample_id": _text,
    "contig_id": _text,
    "truth": _text,
    "contig_length": _nonnegative_int,
}


RowT = TypeVar("RowT")


def read_tsv(path: str | Path, row_type: type[RowT]) -> list[RowT]:
    """Read a TSV with an exact header and report row/column validation errors."""

    path = Path(path)
    parsers = getattr(row_type, "PARSERS")
    expected = list(parsers)
    rows: list[RowT] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != expected:
            raise SchemaError(
                f"{path}: expected header {expected}, found {reader.fieldnames}"
            )
        for line_number, raw in enumerate(reader, start=2):
            parsed: dict[str, object] = {}
            for name, parser in parsers.items():
                try:
                    parsed[name] = parser(raw[name])
                except (TypeError, ValueError) as exc:
                    raise SchemaError(
                        f"{path}:{line_number}: invalid {name}: {exc}"
                    ) from exc
            try:
                rows.append(row_type(**parsed))
            except (TypeError, ValueError) as exc:
                raise SchemaError(f"{path}:{line_number}: {exc}") from exc
    if not rows:
        raise SchemaError(f"{path}: table must contain at least one data row")
    return rows


def assert_unique(rows: Iterable[object], key_fields: tuple[str, ...]) -> None:
    """Reject duplicate natural keys before joins can inflate denominators."""

    seen: set[tuple[object, ...]] = set()
    for row in rows:
        key = tuple(getattr(row, field) for field in key_fields)
        if key in seen:
            raise SchemaError(f"duplicate key {key_fields}={key}")
        seen.add(key)


def row_field_names(row_type: type[object]) -> list[str]:
    """Expose serializable fields without the class-level parser registry."""

    return [field.name for field in fields(row_type) if field.name != "PARSERS"]
