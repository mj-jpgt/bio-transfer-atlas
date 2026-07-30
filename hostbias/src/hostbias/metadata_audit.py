"""Generate commit-safe evidence that live ENA metadata reproduces the freeze."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path
from typing import Callable, Sequence

from .data_manifest import (
    ENA_FIELDS,
    MANIFEST_FIELDS,
    ManifestError,
    canonical_tsv,
    eligible_rows,
    fetch_ena,
    read_tsv,
    select_arm,
    sha256_bytes,
    validate_manifest,
    write_snapshot,
)

Fetcher = Callable[[str], list[dict[str, str]]]


def audit_live_metadata(
    frozen_manifest_path: Path,
    arm_projects: Sequence[tuple[str, str]],
    *,
    primary: int = 20,
    reserves: int = 10,
    fetcher: Fetcher = fetch_ena,
    snapshot_dir: Path | None = None,
    checked_at_utc: str | None = None,
) -> dict[str, object]:
    """Fetch ENA, rebuild selection, and return only aggregate/hash evidence."""
    if not arm_projects:
        raise ManifestError("at least one ARM=PROJECT specification is required")
    arm_names = [arm for arm, _ in arm_projects]
    if len(arm_names) != len(set(arm_names)):
        raise ManifestError("arm specifications must be unique")

    frozen = read_tsv(frozen_manifest_path)
    validate_manifest(frozen, primary=primary, reserves=reserves)
    live_selected: list[dict[str, object]] = []
    arm_evidence: dict[str, object] = {}
    for arm, project in arm_projects:
        rows = fetcher(project)
        snapshot_payload = canonical_tsv(
            sorted(rows, key=lambda row: row["run_accession"]), ENA_FIELDS
        )
        if snapshot_dir is not None:
            write_snapshot(
                snapshot_dir / f"{project}.ena.tsv",
                sorted(rows, key=lambda row: row["run_accession"]),
                ENA_FIELDS,
            )
        eligible = eligible_rows(rows)
        selected = select_arm(
            rows,
            arm=arm,
            project=project,
            primary=primary,
            reserves=reserves,
        )
        live_selected.extend(selected)
        frozen_arm = [row for row in frozen if row["arm"] == arm]
        live_arm_strings = [
            {field: str(row[field]) for field in MANIFEST_FIELDS} for row in selected
        ]
        arm_evidence[arm] = {
            "bioproject": project,
            "response_rows": len(rows),
            "eligible_unique_rows": len(eligible),
            "canonical_snapshot_sha256": sha256_bytes(snapshot_payload),
            "selection_matches_frozen": live_arm_strings == frozen_arm,
            "selected_primary": primary,
            "selected_reserve": reserves,
        }

    live_strings = [
        {field: str(row[field]) for field in MANIFEST_FIELDS} for row in live_selected
    ]
    validate_manifest(live_strings, primary=primary, reserves=reserves)
    frozen_payload = canonical_tsv(frozen, MANIFEST_FIELDS)
    live_payload = canonical_tsv(live_selected, MANIFEST_FIELDS)
    exact_match = live_payload == frozen_payload
    if checked_at_utc is None:
        checked_at_utc = (
            dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
        )
    return {
        "schema_version": 1,
        "checked_at_utc": checked_at_utc,
        "valid": exact_match
        and all(
            bool(evidence["selection_matches_frozen"])
            for evidence in arm_evidence.values()
        ),
        "frozen_manifest_sha256": sha256_bytes(frozen_payload),
        "live_manifest_sha256": sha256_bytes(live_payload),
        "manifest_exact_match": exact_match,
        "arms": arm_evidence,
    }


def _arm_project(value: str) -> tuple[str, str]:
    parts = value.split("=", 1)
    if len(parts) != 2 or not all(parts):
        raise argparse.ArgumentTypeError("expected ARM=PROJECT")
    return parts[0], parts[1]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument(
        "--arm",
        required=True,
        action="append",
        type=_arm_project,
        metavar="ARM=PROJECT",
    )
    parser.add_argument("--primary", type=int, default=20)
    parser.add_argument("--reserves", type=int, default=10)
    parser.add_argument("--snapshot-dir", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    report = audit_live_metadata(
        args.manifest,
        args.arm,
        primary=args.primary,
        reserves=args.reserves,
        snapshot_dir=args.snapshot_dir,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps({"valid": report["valid"], "output": str(args.output)}))
    return 0 if report["valid"] else 5


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ManifestError as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
