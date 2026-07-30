"""Reproduce and validate the frozen HPRC donor selection from public metadata."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .data_manifest import ManifestError, read_tsv

SUPERPOPULATIONS = ("AFR", "AMR", "EAS", "EUR", "SAS")
DONOR_FIELDS = (
    "donor_id",
    "biosample_accession",
    "population_code",
    "superpopulation",
    "panel_role",
    "selection_rank",
    "haplotype_count",
    "metadata_basis",
)
PANEL_FIELDS = (
    "reference_id",
    "kind",
    "donor_id",
    "haplotype",
    "population_code",
    "superpopulation",
    "panel_role",
    "genbank_accession",
    "source_url",
    "source_md5_url",
    "expected_bytes",
    "expected_md5",
    "expected_sha256",
    "include_in_union",
    "local_filename",
)


def read_csv(path: Path, *, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    if not rows:
        raise ManifestError(f"{path} contains no metadata rows")
    return rows


def _require_fields(
    rows: Sequence[Mapping[str, str]], fields: Iterable[str], context: str
) -> None:
    missing = set(fields).difference(rows[0] if rows else ())
    if missing:
        raise ManifestError(f"{context} missing fields: {sorted(missing)}")


def _https_hprc_url(value: str) -> str:
    prefix = "s3://human-pangenomics/"
    if value.startswith(prefix):
        return value.replace(
            prefix,
            "https://s3-us-west-2.amazonaws.com/human-pangenomics/",
            1,
        )
    if value.startswith("https://"):
        return value
    raise ManifestError(f"unsupported HPRC assembly URL: {value!r}")


def build_balanced_donor_catalog(
    assemblies: Sequence[Mapping[str, str]],
    samples: Sequence[Mapping[str, str]],
    igsr_panel: Sequence[Mapping[str, str]],
    *,
    primary_per_group: int = 1,
    holdout_per_group: int = 1,
) -> list[dict[str, str]]:
    """Select donors only after exact sample and population-code joins."""
    _require_fields(
        assemblies,
        (
            "sample_id",
            "haplotype",
            "genbank_accession",
            "assembly_md5",
            "assembly",
        ),
        "HPRC assembly metadata",
    )
    _require_fields(
        samples,
        ("sample_id", "biosample_id", "population_abbreviation"),
        "HPRC sample metadata",
    )
    _require_fields(
        igsr_panel,
        ("sample", "pop", "super_pop"),
        "IGSR panel metadata",
    )
    if primary_per_group <= 0 or holdout_per_group <= 0:
        raise ManifestError("primary and holdout counts must be positive")

    sample_by_id = {row["sample_id"]: row for row in samples}
    if len(sample_by_id) != len(samples):
        raise ManifestError("HPRC sample metadata has duplicate sample_id values")
    igsr_by_id = {row["sample"]: row for row in igsr_panel}
    if len(igsr_by_id) != len(igsr_panel):
        raise ManifestError("IGSR panel metadata has duplicate sample values")
    assemblies_by_donor: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in assemblies:
        assemblies_by_donor[row["sample_id"]].append(row)

    eligible: dict[str, list[tuple[str, Mapping[str, str]]]] = defaultdict(list)
    for donor, donor_assemblies in assemblies_by_donor.items():
        sample = sample_by_id.get(donor)
        igsr = igsr_by_id.get(donor)
        if sample is None or igsr is None:
            continue
        if sample["population_abbreviation"] != igsr["pop"]:
            continue
        if igsr["super_pop"] not in SUPERPOPULATIONS:
            continue
        if len(donor_assemblies) != 2:
            continue
        if {row["haplotype"] for row in donor_assemblies} != {"1", "2"}:
            continue
        if any(
            not row["genbank_accession"]
            or not row["assembly"]
            or not row["assembly_md5"]
            for row in donor_assemblies
        ):
            continue
        if not sample["biosample_id"]:
            continue
        eligible[igsr["super_pop"]].append((donor, sample))

    result: list[dict[str, str]] = []
    required = primary_per_group + holdout_per_group
    for superpopulation in SUPERPOPULATIONS:
        candidates = sorted(eligible[superpopulation], key=lambda item: item[0])
        if len(candidates) < required:
            raise ManifestError(
                f"{superpopulation}: need {required} eligible HPRC donors; "
                f"found {len(candidates)}"
            )
        for rank, (donor, sample) in enumerate(candidates[:required], start=1):
            result.append(
                {
                    "donor_id": donor,
                    "biosample_accession": sample["biosample_id"],
                    "population_code": sample["population_abbreviation"],
                    "superpopulation": superpopulation,
                    "panel_role": (
                        "primary" if rank <= primary_per_group else "holdout"
                    ),
                    "selection_rank": str(rank),
                    "haplotype_count": "2",
                    "metadata_basis": "HPRC+IGSR exact donor join",
                }
            )
    validate_donor_catalog(
        result,
        primary_per_group=primary_per_group,
        holdout_per_group=holdout_per_group,
    )
    return result


def validate_donor_catalog(
    rows: Sequence[Mapping[str, str]],
    *,
    primary_per_group: int = 1,
    holdout_per_group: int = 1,
) -> dict[str, object]:
    _require_fields(rows, DONOR_FIELDS, "donor catalog")
    donors = [row["donor_id"] for row in rows]
    if len(donors) != len(set(donors)):
        raise ManifestError("donor catalog contains duplicate donors")
    if {row["superpopulation"] for row in rows} != set(SUPERPOPULATIONS):
        raise ManifestError("donor catalog must contain AFR, AMR, EAS, EUR, and SAS")
    for superpopulation in SUPERPOPULATIONS:
        group = [row for row in rows if row["superpopulation"] == superpopulation]
        roles = Counter(row["panel_role"] for row in group)
        if roles != {
            "primary": primary_per_group,
            "holdout": holdout_per_group,
        }:
            raise ManifestError(f"{superpopulation}: unbalanced donor roles {dict(roles)}")
        expected_ranks = list(range(1, primary_per_group + holdout_per_group + 1))
        observed_ranks = [int(row["selection_rank"]) for row in group]
        if observed_ranks != expected_ranks:
            raise ManifestError(f"{superpopulation}: invalid selection ranks")
        if any(row["haplotype_count"] != "2" for row in group):
            raise ManifestError(f"{superpopulation}: donors require two haplotypes")
    return {
        "valid": True,
        "donors": len(rows),
        "superpopulations": list(SUPERPOPULATIONS),
        "primary_per_group": primary_per_group,
        "holdout_per_group": holdout_per_group,
    }


def validate_panel_manifest(
    panel: Sequence[Mapping[str, str]],
    donors: Sequence[Mapping[str, str]],
) -> dict[str, object]:
    _require_fields(panel, PANEL_FIELDS, "competitive panel manifest")
    validate_donor_catalog(donors)
    donor_by_id = {row["donor_id"]: row for row in donors}
    reference_ids = [row["reference_id"] for row in panel]
    filenames = [row["local_filename"] for row in panel]
    if len(reference_ids) != len(set(reference_ids)):
        raise ManifestError("competitive panel contains duplicate reference_id values")
    if len(filenames) != len(set(filenames)):
        raise ManifestError("competitive panel contains duplicate filenames")
    chm13 = [row for row in panel if row["kind"] == "chm13"]
    if len(chm13) != 1 or chm13[0]["include_in_union"] != "true":
        raise ManifestError("competitive panel requires one union-included CHM13")
    hprc = [row for row in panel if row["kind"] == "hprc_assembly"]
    if len(hprc) != len(donors) * 2:
        raise ManifestError("competitive panel requires two assemblies per donor")
    for row in hprc:
        donor = donor_by_id.get(row["donor_id"])
        if donor is None:
            raise ManifestError(f"{row['reference_id']}: donor is not frozen")
        for field in ("population_code", "superpopulation", "panel_role"):
            if row[field] != donor[field]:
                raise ManifestError(
                    f"{row['reference_id']}: {field} disagrees with donor catalog"
                )
        if row["haplotype"] not in {"1", "2"}:
            raise ManifestError(f"{row['reference_id']}: invalid haplotype")
        should_include = donor["panel_role"] == "primary"
        if (row["include_in_union"] == "true") != should_include:
            raise ManifestError(
                f"{row['reference_id']}: union membership disagrees with donor role"
            )
        try:
            expected_bytes = int(row["expected_bytes"])
        except ValueError as exc:
            raise ManifestError(
                f"{row['reference_id']}: expected_bytes must be an integer"
            ) from exc
        if expected_bytes <= 0:
            raise ManifestError(f"{row['reference_id']}: expected_bytes must be positive")
        if len(row["expected_md5"]) != 32:
            raise ManifestError(f"{row['reference_id']}: invalid expected MD5")
    primary_counts = Counter(
        row["superpopulation"]
        for row in hprc
        if row["include_in_union"] == "true"
    )
    if primary_counts != Counter({group: 2 for group in SUPERPOPULATIONS}):
        raise ManifestError("union HPRC haplotypes are not balanced")
    return {
        "valid": True,
        "references": len(panel),
        "hprc_donors": len(donors),
        "union_references": sum(row["include_in_union"] == "true" for row in panel),
        "expected_download_bytes": sum(int(row["expected_bytes"]) for row in panel),
        "expected_union_input_bytes": sum(
            int(row["expected_bytes"])
            for row in panel
            if row["include_in_union"] == "true"
        ),
    }


def verify_frozen_selection(
    *,
    assembly_metadata: Path,
    sample_metadata: Path,
    igsr_panel: Path,
    frozen_donors: Path,
    panel_manifest: Path,
) -> dict[str, object]:
    reproduced = build_balanced_donor_catalog(
        read_csv(assembly_metadata),
        read_csv(sample_metadata),
        read_csv(igsr_panel, delimiter="\t"),
    )
    frozen = read_tsv(frozen_donors)
    if reproduced != frozen:
        raise ManifestError(
            "authoritative metadata no longer reproduces frozen donor selection"
        )
    panel_report = validate_panel_manifest(read_tsv(panel_manifest), frozen)
    return {
        "valid": True,
        "selection_exact_match": True,
        "donor_catalog": validate_donor_catalog(frozen),
        "panel_manifest": panel_report,
    }


def validate_panel_against_assemblies(
    panel: Sequence[Mapping[str, str]],
    assemblies: Sequence[Mapping[str, str]],
) -> None:
    """Verify each frozen HPRC row against its exact authoritative assembly row."""
    by_key = {
        (row["sample_id"], row["haplotype"]): row
        for row in assemblies
        if row.get("sample_id") and row.get("haplotype")
    }
    for frozen in panel:
        if frozen["kind"] != "hprc_assembly":
            continue
        key = (frozen["donor_id"], frozen["haplotype"])
        source = by_key.get(key)
        if source is None:
            raise ManifestError(f"{frozen['reference_id']}: missing assembly metadata")
        expected_values = {
            "genbank_accession": source["genbank_accession"],
            "source_url": _https_hprc_url(source["assembly"]),
            "source_md5_url": _https_hprc_url(source["assembly_md5"]),
        }
        for field, expected in expected_values.items():
            if frozen[field] != expected:
                raise ManifestError(
                    f"{frozen['reference_id']}: {field} disagrees with HPRC metadata"
                )

