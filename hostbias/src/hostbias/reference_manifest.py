"""Build and verify ancestry-balanced human reference manifests."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.parse
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .data_manifest import ManifestError, read_tsv, write_snapshot

CATALOG_FIELDS = (
    "donor_id",
    "population_group",
    "assembly_url",
    "sha256",
)
REFERENCE_FIELDS = (
    "reference_id",
    "kind",
    "population_group",
    "donor_id",
    "source_url",
    "expected_sha256",
    "included_in_primary",
    "held_out_control",
    "local_filename",
)
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def _validate_sha256(value: str, context: str) -> str:
    normalized = value.strip().lower()
    if not SHA256_PATTERN.fullmatch(normalized):
        raise ManifestError(f"{context}: sha256 must be 64 hexadecimal characters")
    return normalized


def _filename_from_url(url: str, context: str) -> str:
    filename = Path(urllib.parse.urlparse(url).path).name
    if not filename or filename in {".", ".."}:
        raise ManifestError(f"{context}: source URL has no filename")
    return filename


def build_balanced_panel(
    catalog: Sequence[Mapping[str, str]],
    *,
    per_group: int,
    holdout_per_group: int = 1,
    chm13_url: str,
    chm13_sha256: str,
) -> list[dict[str, str]]:
    """Select equal primary/held-out donor counts per declared population group."""
    if per_group <= 0 or holdout_per_group <= 0:
        raise ManifestError("per-group and holdout-per-group must be positive")
    if not catalog:
        raise ManifestError("HPRC catalog is empty")
    missing = set(CATALOG_FIELDS).difference(catalog[0])
    if missing:
        raise ManifestError(f"HPRC catalog missing fields: {sorted(missing)}")

    donors = [row["donor_id"] for row in catalog]
    if any(not donor for donor in donors):
        raise ManifestError("HPRC catalog contains an empty donor_id")
    if len(donors) != len(set(donors)):
        raise ManifestError("HPRC catalog donor_id values are not unique")

    panel: list[dict[str, str]] = [
        {
            "reference_id": "chm13v2.0",
            "kind": "chm13",
            "population_group": "reference",
            "donor_id": "CHM13",
            "source_url": chm13_url,
            "expected_sha256": _validate_sha256(chm13_sha256, "CHM13"),
            "included_in_primary": "true",
            "held_out_control": "false",
            "local_filename": _filename_from_url(chm13_url, "CHM13"),
        }
    ]
    groups = sorted({row["population_group"] for row in catalog if row["population_group"]})
    if not groups:
        raise ManifestError("HPRC catalog has no population groups")
    required = per_group + holdout_per_group
    for group in groups:
        candidates = sorted(
            (row for row in catalog if row["population_group"] == group),
            key=lambda row: row["donor_id"],
        )
        if len(candidates) < required:
            raise ManifestError(
                f"{group}: need {required} donors for primary plus holdout; "
                f"found {len(candidates)}"
            )
        for position, row in enumerate(candidates[:required], start=1):
            held_out = position > per_group
            donor = row["donor_id"]
            url = row["assembly_url"]
            panel.append(
                {
                    "reference_id": f"hprc_{donor}",
                    "kind": "hprc_assembly",
                    "population_group": group,
                    "donor_id": donor,
                    "source_url": url,
                    "expected_sha256": _validate_sha256(
                        row["sha256"], f"HPRC donor {donor}"
                    ),
                    "included_in_primary": str(not held_out).lower(),
                    "held_out_control": str(held_out).lower(),
                    "local_filename": _filename_from_url(url, f"HPRC donor {donor}"),
                }
            )
    validate_reference_manifest(
        panel, per_group=per_group, holdout_per_group=holdout_per_group
    )
    return panel


def validate_reference_manifest(
    rows: Sequence[Mapping[str, str]],
    *,
    per_group: int | None = None,
    holdout_per_group: int | None = None,
) -> dict[str, object]:
    if not rows:
        raise ManifestError("reference manifest is empty")
    missing = set(REFERENCE_FIELDS).difference(rows[0])
    if missing:
        raise ManifestError(f"reference manifest missing fields: {sorted(missing)}")
    reference_ids = [row["reference_id"] for row in rows]
    if len(reference_ids) != len(set(reference_ids)):
        raise ManifestError("reference_id values are not unique")
    filenames = [row["local_filename"] for row in rows]
    if len(filenames) != len(set(filenames)):
        raise ManifestError("local_filename values are not unique")
    for row in rows:
        _validate_sha256(row["expected_sha256"], row["reference_id"])
        if row["included_in_primary"] not in {"true", "false"}:
            raise ManifestError(
                f"{row['reference_id']}: included_in_primary must be true or false"
            )
        if row["held_out_control"] not in {"true", "false"}:
            raise ManifestError(
                f"{row['reference_id']}: held_out_control must be true or false"
            )
        if (
            row["included_in_primary"] == "true"
            and row["held_out_control"] == "true"
        ):
            raise ManifestError(
                f"{row['reference_id']}: a held-out control cannot enter primary panel"
            )

    chm13 = [row for row in rows if row["kind"] == "chm13"]
    if len(chm13) != 1 or chm13[0]["included_in_primary"] != "true":
        raise ManifestError("manifest must contain one primary CHM13 reference")
    hprc = [row for row in rows if row["kind"] == "hprc_assembly"]
    primary_counts = Counter(
        row["population_group"]
        for row in hprc
        if row["included_in_primary"] == "true"
    )
    holdout_counts = Counter(
        row["population_group"]
        for row in hprc
        if row["held_out_control"] == "true"
    )
    if not primary_counts or len(set(primary_counts.values())) != 1:
        raise ManifestError("HPRC primary panel is not balanced across groups")
    if set(primary_counts) != set(holdout_counts):
        raise ManifestError("every primary population group requires held-out donors")
    if len(set(holdout_counts.values())) != 1:
        raise ManifestError("HPRC held-out panel is not balanced across groups")
    if per_group is not None and set(primary_counts.values()) != {per_group}:
        raise ManifestError("primary donor count does not match per_group")
    if (
        holdout_per_group is not None
        and set(holdout_counts.values()) != {holdout_per_group}
    ):
        raise ManifestError("held-out donor count does not match holdout_per_group")
    return {
        "valid": True,
        "references": len(rows),
        "population_groups": sorted(primary_counts),
        "primary_hprc_per_group": next(iter(primary_counts.values())),
        "held_out_hprc_per_group": next(iter(holdout_counts.values())),
    }


def sha256_file(path: Path, chunk_bytes: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_bytes):
            digest.update(chunk)
    return digest.hexdigest()


def verify_downloads(
    rows: Iterable[Mapping[str, str]], reference_root: Path
) -> dict[str, object]:
    checks = []
    for row in rows:
        path = reference_root / row["local_filename"]
        expected = _validate_sha256(row["expected_sha256"], row["reference_id"])
        if path.is_file():
            observed = sha256_file(path)
            status = "ok" if observed == expected else "checksum_mismatch"
        else:
            observed = None
            status = "missing"
        checks.append(
            {
                "reference_id": row["reference_id"],
                "path": str(path),
                "expected_sha256": expected,
                "observed_sha256": observed,
                "status": status,
            }
        )
    return {
        "valid": all(check["status"] == "ok" for check in checks),
        "checks": checks,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    build = commands.add_parser("build", help="select and freeze a balanced panel")
    build.add_argument("--hprc-catalog", required=True, type=Path)
    build.add_argument("--per-group", required=True, type=int)
    build.add_argument("--holdout-per-group", default=1, type=int)
    build.add_argument("--chm13-url", required=True)
    build.add_argument("--chm13-sha256", required=True)
    build.add_argument("--output", required=True, type=Path)

    validate = commands.add_parser("validate", help="validate panel balance")
    validate.add_argument("--manifest", required=True, type=Path)

    verify = commands.add_parser("verify", help="verify downloaded reference checksums")
    verify.add_argument("--manifest", required=True, type=Path)
    verify.add_argument("--reference-root", required=True, type=Path)
    verify.add_argument("--report", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "build":
        rows = build_balanced_panel(
            read_tsv(args.hprc_catalog),
            per_group=args.per_group,
            holdout_per_group=args.holdout_per_group,
            chm13_url=args.chm13_url,
            chm13_sha256=args.chm13_sha256,
        )
        digest = write_snapshot(args.output, rows, REFERENCE_FIELDS)
        print(json.dumps({"references": len(rows), "sha256": digest}))
        return 0
    if args.command == "validate":
        report = validate_reference_manifest(read_tsv(args.manifest))
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0

    report = verify_downloads(read_tsv(args.manifest), args.reference_root)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps({"valid": report["valid"], "report": str(args.report)}))
    return 0 if report["valid"] else 4


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ManifestError as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)

