"""Evaluate the preregistered Stage 0 sentinel eligibility check."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Mapping, Sequence

from .data_manifest import ManifestError, read_tsv

METRIC_FIELDS = (
    "run_accession",
    "arm",
    "streamed_spots",
    "grch38_mapped_pairs",
    "metadata_ok",
    "checksum_ok",
)
TRUE_VALUES = {"1", "true", "yes"}


def _as_bool(value: str, *, field: str, run: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in TRUE_VALUES | {"0", "false", "no"}:
        raise ManifestError(f"{run}: {field} must be true or false")
    return normalized in TRUE_VALUES


def _as_positive_int(value: str, *, field: str, run: str) -> int:
    try:
        result = int(value)
    except ValueError as exc:
        raise ManifestError(f"{run}: {field} must be an integer") from exc
    if result <= 0:
        raise ManifestError(f"{run}: {field} must be positive")
    return result


def evaluate_sentinels(
    manifest: Sequence[Mapping[str, str]],
    metrics: Sequence[Mapping[str, str]],
    *,
    sentinels_per_arm: int = 3,
    minimum_mapped_pairs_per_million: float = 100.0,
    minimum_passing_runs: int = 2,
) -> dict[str, object]:
    if not metrics:
        raise ManifestError("sentinel metrics are empty")
    missing_fields = set(METRIC_FIELDS).difference(metrics[0])
    if missing_fields:
        raise ManifestError(f"sentinel metrics missing fields: {sorted(missing_fields)}")
    by_accession = {row["run_accession"]: row for row in metrics}
    if len(by_accession) != len(metrics):
        raise ManifestError("sentinel metrics contain duplicate run_accession values")

    arm_reports: dict[str, object] = {}
    for arm in sorted({row["arm"] for row in manifest}):
        expected = [
            row
            for row in manifest
            if row["arm"] == arm
            and row["role"] == "primary"
            and int(row["rank"]) <= sentinels_per_arm
        ]
        if len(expected) != sentinels_per_arm:
            raise ManifestError(
                f"{arm}: expected {sentinels_per_arm} frozen sentinel runs, "
                f"found {len(expected)}"
            )
        run_reports = []
        for frozen in expected:
            run = frozen["run_accession"]
            if run not in by_accession:
                raise ManifestError(f"{arm}: missing sentinel metrics for {run}")
            observed = by_accession[run]
            if observed["arm"] != arm:
                raise ManifestError(f"{run}: arm does not match frozen manifest")
            spots = _as_positive_int(
                observed["streamed_spots"], field="streamed_spots", run=run
            )
            try:
                mapped = int(observed["grch38_mapped_pairs"])
            except ValueError as exc:
                raise ManifestError(
                    f"{run}: grch38_mapped_pairs must be an integer"
                ) from exc
            if mapped < 0 or mapped > spots:
                raise ManifestError(
                    f"{run}: grch38_mapped_pairs must be between 0 and streamed_spots"
                )
            mapped_per_million = mapped / spots * 1_000_000
            metadata_ok = _as_bool(observed["metadata_ok"], field="metadata_ok", run=run)
            checksum_ok = _as_bool(observed["checksum_ok"], field="checksum_ok", run=run)
            host_signal_pass = (
                mapped_per_million >= minimum_mapped_pairs_per_million
            )
            run_reports.append(
                {
                    "run_accession": run,
                    "mapped_pairs_per_million": mapped_per_million,
                    "host_signal_pass": host_signal_pass,
                    "metadata_ok": metadata_ok,
                    "checksum_ok": checksum_ok,
                }
            )
        passing_host_runs = sum(
            bool(run["host_signal_pass"]) for run in run_reports
        )
        all_provenance_ok = all(
            bool(run["metadata_ok"]) and bool(run["checksum_ok"])
            for run in run_reports
        )
        arm_reports[arm] = {
            "eligible": (
                passing_host_runs >= minimum_passing_runs and all_provenance_ok
            ),
            "passing_host_runs": passing_host_runs,
            "required_passing_host_runs": minimum_passing_runs,
            "all_metadata_and_checksums_ok": all_provenance_ok,
            "runs": run_reports,
        }
    return {
        "eligible": all(bool(report["eligible"]) for report in arm_reports.values()),
        "rule": {
            "sentinels_per_arm": sentinels_per_arm,
            "minimum_mapped_pairs_per_million": minimum_mapped_pairs_per_million,
            "minimum_passing_runs": minimum_passing_runs,
            "require_all_metadata_and_checksums_ok": True,
        },
        "arms": arm_reports,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sentinels-per-arm", type=int, default=3)
    parser.add_argument("--minimum-mapped-pairs-per-million", type=float, default=100)
    parser.add_argument("--minimum-passing-runs", type=int, default=2)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    report = evaluate_sentinels(
        read_tsv(args.manifest),
        read_tsv(args.metrics),
        sentinels_per_arm=args.sentinels_per_arm,
        minimum_mapped_pairs_per_million=args.minimum_mapped_pairs_per_million,
        minimum_passing_runs=args.minimum_passing_runs,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps({"eligible": report["eligible"], "output": str(args.output)}))
    return 0 if report["eligible"] else 3


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ManifestError as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)

