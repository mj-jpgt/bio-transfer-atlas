"""Deterministic ENA metadata freezing and Stage 0 sample selection."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping, Sequence

ENA_ENDPOINT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = (
    "run_accession",
    "study_accession",
    "sample_accession",
    "library_layout",
    "instrument_platform",
    "library_strategy",
    "base_count",
    "read_count",
    "fastq_ftp",
    "fastq_md5",
)
RUNTIME_ENA_FIELDS = ENA_FIELDS + ("fastq_bytes",)
MANIFEST_FIELDS = (
    "run_accession",
    "arm",
    "bioproject",
    "sample_accession",
    "role",
    "rank",
    "library_layout",
    "instrument_platform",
    "library_strategy",
    "base_count",
    "read_count",
)


class ManifestError(ValueError):
    """Raised when metadata or a frozen manifest violates its contract."""


def _read_tsv_text(text: str) -> list[dict[str, str]]:
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        raise ManifestError("TSV contains no data rows")
    return rows


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ManifestError(f"{path} contains no data rows")
    return rows


def canonical_tsv(rows: Iterable[Mapping[str, object]], fields: Sequence[str]) -> bytes:
    """Return stable LF-delimited TSV bytes with a fixed column order."""
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in fields})
    return buffer.getvalue().encode("utf-8")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def write_snapshot(
    output: Path,
    rows: Iterable[Mapping[str, object]],
    fields: Sequence[str] = ENA_FIELDS,
) -> str:
    payload = canonical_tsv(rows, fields)
    digest = sha256_bytes(payload)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(payload)
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{digest}  {output.name}\n", encoding="ascii", newline="\n"
    )
    return digest


def fetch_ena(
    project: str,
    timeout_seconds: int = 120,
    fields: Sequence[str] = ENA_FIELDS,
) -> list[dict[str, str]]:
    params = urllib.parse.urlencode(
        {
            "accession": project,
            "result": "read_run",
            "fields": ",".join(fields),
            "format": "tsv",
            "download": "false",
        }
    )
    request = urllib.request.Request(
        f"{ENA_ENDPOINT}?{params}",
        headers={"User-Agent": "hostbias-metadata/1"},
    )
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
        text = response.read().decode("utf-8")
    rows = _read_tsv_text(text)
    missing = set(fields).difference(rows[0])
    if missing:
        raise ManifestError(f"ENA response missing fields: {sorted(missing)}")
    return sorted(rows, key=lambda row: row["run_accession"])


def _positive_int(row: Mapping[str, str], field: str) -> int:
    try:
        value = int(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise ManifestError(
            f"{row.get('run_accession', '<unknown>')}: invalid {field}"
        ) from exc
    if value < 0:
        raise ManifestError(
            f"{row.get('run_accession', '<unknown>')}: negative {field}"
        )
    return value


def eligible_rows(rows: Iterable[Mapping[str, str]]) -> list[Mapping[str, str]]:
    """Apply the frozen technical filters and one-run-per-BioSample rule."""
    technical = [
        row
        for row in rows
        if row.get("library_layout") == "PAIRED"
        and row.get("instrument_platform") == "ILLUMINA"
        and row.get("library_strategy") == "WGS"
    ]
    ranked = sorted(
        technical,
        key=lambda row: (-_positive_int(row, "base_count"), row["run_accession"]),
    )
    seen_samples: set[str] = set()
    unique: list[Mapping[str, str]] = []
    for row in ranked:
        sample = row.get("sample_accession", "")
        if not sample:
            raise ManifestError(f"{row.get('run_accession')}: missing sample_accession")
        if sample not in seen_samples:
            seen_samples.add(sample)
            unique.append(row)
    return unique


def select_arm(
    rows: Iterable[Mapping[str, str]],
    *,
    arm: str,
    project: str,
    primary: int = 20,
    reserves: int = 10,
) -> list[dict[str, object]]:
    selected = eligible_rows(rows)[: primary + reserves]
    if len(selected) < primary + reserves:
        raise ManifestError(
            f"{arm}: need {primary + reserves} eligible unique samples; "
            f"found {len(selected)}"
        )
    result: list[dict[str, object]] = []
    for position, row in enumerate(selected, start=1):
        result.append(
            {
                "run_accession": row["run_accession"],
                "arm": arm,
                "bioproject": project,
                "sample_accession": row["sample_accession"],
                "role": "primary" if position <= primary else "reserve",
                "rank": position,
                "library_layout": row["library_layout"],
                "instrument_platform": row["instrument_platform"],
                "library_strategy": row["library_strategy"],
                "base_count": _positive_int(row, "base_count"),
                "read_count": _positive_int(row, "read_count"),
            }
        )
    return result


def validate_manifest(
    rows: Sequence[Mapping[str, str]],
    *,
    primary: int = 20,
    reserves: int = 10,
) -> dict[str, object]:
    missing = set(MANIFEST_FIELDS).difference(rows[0] if rows else ())
    if missing:
        raise ManifestError(f"manifest missing fields: {sorted(missing)}")
    runs = [row["run_accession"] for row in rows]
    if len(runs) != len(set(runs)):
        raise ManifestError("run_accession values are not unique")

    arms: dict[str, dict[str, object]] = {}
    for arm in sorted({row["arm"] for row in rows}):
        arm_rows = [row for row in rows if row["arm"] == arm]
        expected_count = primary + reserves
        if len(arm_rows) != expected_count:
            raise ManifestError(
                f"{arm}: expected {expected_count} rows, found {len(arm_rows)}"
            )
        expected_ranks = list(range(1, expected_count + 1))
        ranks = [_positive_int(row, "rank") for row in arm_rows]
        if ranks != expected_ranks:
            raise ManifestError(f"{arm}: ranks must be ordered 1..{expected_count}")
        roles = Counter(row["role"] for row in arm_rows)
        if roles != {"primary": primary, "reserve": reserves}:
            raise ManifestError(f"{arm}: invalid role counts {dict(roles)}")
        if any(row["role"] != ("primary" if i <= primary else "reserve")
               for i, row in enumerate(arm_rows, start=1)):
            raise ManifestError(f"{arm}: role does not match rank")
        bases = [_positive_int(row, "base_count") for row in arm_rows]
        accessions = [row["run_accession"] for row in arm_rows]
        observed_keys = list(zip((-value for value in bases), accessions))
        if observed_keys != sorted(observed_keys):
            raise ManifestError(
                f"{arm}: rows are not ranked by base_count desc/accession asc"
            )
        samples = [row["sample_accession"] for row in arm_rows]
        if len(samples) != len(set(samples)):
            raise ManifestError(f"{arm}: sample_accession values are not unique")
        arms[arm] = {
            "rows": len(arm_rows),
            "primary": roles["primary"],
            "reserve": roles["reserve"],
            "bioproject": sorted({row["bioproject"] for row in arm_rows}),
        }
    return {"valid": True, "total_rows": len(rows), "arms": arms}


def _parse_arm_spec(spec: str) -> tuple[str, str, Path]:
    parts = spec.split("=", 2)
    if len(parts) != 3 or not all(parts):
        raise argparse.ArgumentTypeError("expected ARM=PROJECT=METADATA.tsv")
    return parts[0], parts[1], Path(parts[2])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    fetch = commands.add_parser("fetch-ena", help="freeze a canonical ENA snapshot")
    fetch.add_argument("--project", required=True)
    fetch.add_argument("--output", required=True, type=Path)

    fetch_runtime = commands.add_parser(
        "fetch-runtime-ena",
        help="freeze a canonical ENA snapshot with FASTQ sizes",
    )
    fetch_runtime.add_argument("--project", required=True)
    fetch_runtime.add_argument("--output", required=True, type=Path)

    select = commands.add_parser("select", help="select primary and reserve runs")
    select.add_argument(
        "--arm",
        action="append",
        required=True,
        type=_parse_arm_spec,
        metavar="ARM=PROJECT=METADATA.tsv",
    )
    select.add_argument("--primary", type=int, default=20)
    select.add_argument("--reserves", type=int, default=10)
    select.add_argument("--output", required=True, type=Path)

    validate = commands.add_parser("validate", help="validate a frozen manifest")
    validate.add_argument("--manifest", required=True, type=Path)
    validate.add_argument("--primary", type=int, default=20)
    validate.add_argument("--reserves", type=int, default=10)
    validate.add_argument("--report", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command in {"fetch-ena", "fetch-runtime-ena"}:
        fields = RUNTIME_ENA_FIELDS if args.command == "fetch-runtime-ena" else ENA_FIELDS
        rows = fetch_ena(args.project, fields=fields)
        digest = write_snapshot(args.output, rows, fields)
        print(json.dumps({"project": args.project, "rows": len(rows), "sha256": digest}))
        return 0
    if args.command == "select":
        selected: list[dict[str, object]] = []
        for arm, project, path in args.arm:
            selected.extend(
                select_arm(
                    read_tsv(path),
                    arm=arm,
                    project=project,
                    primary=args.primary,
                    reserves=args.reserves,
                )
            )
        validate_manifest(
            [{key: str(value) for key, value in row.items()} for row in selected],
            primary=args.primary,
            reserves=args.reserves,
        )
        digest = write_snapshot(args.output, selected, MANIFEST_FIELDS)
        print(json.dumps({"rows": len(selected), "sha256": digest}))
        return 0

    report = validate_manifest(
        read_tsv(args.manifest), primary=args.primary, reserves=args.reserves
    )
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(payload, encoding="utf-8", newline="\n")
    else:
        sys.stdout.write(payload)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ManifestError as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)

