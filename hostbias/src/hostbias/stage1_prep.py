"""Freeze and validate the Stage 1 experimental design before any outcome run.

This deliberately prepares only a public donor manifest and an aggregate
checkpoint.  It does not download reads, build a graph, spike sequences, or
calculate a biological result.  Those actions remain conditional on a PASS
Gate A verdict.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from pathlib import Path
from typing import Any

import yaml

from hostbias.data_manifest import ManifestError, canonical_tsv, read_tsv, sha256_bytes
from hostbias.provenance import write_json_atomic


DONOR_FIELDS = (
    "donor_id",
    "superpopulation",
    "population_code",
    "assembly_source",
    "assembly_accession",
    "selection_rank",
)
REQUIRED_GROUPS = ("AFR", "AMR", "EAS", "EUR", "SAS")
REQUIRED_FRACTIONS = (0.001, 0.01, 0.05, 0.1)


class Stage1PreparationError(ManifestError):
    """Raised when a Stage 1 plan is underspecified or confounded."""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_design(path: Path) -> dict[str, Any]:
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as error:
        raise Stage1PreparationError("stage 1 design is not valid YAML") from error
    if not isinstance(payload, dict):
        raise Stage1PreparationError("stage 1 design must be a mapping")
    return payload


def _validate_design(design: dict[str, Any]) -> int:
    required = {
        "gate_a_status",
        "donors_per_superpopulation",
        "superpopulations",
        "backgrounds",
        "spike_fractions",
        "analyses",
        "allowed_assembly_sources",
    }
    missing = sorted(required.difference(design))
    if missing:
        raise Stage1PreparationError(f"stage 1 design missing fields {missing}")
    if design["gate_a_status"] != "PASS_REQUIRED":
        raise Stage1PreparationError("stage 1 design must require a Gate A PASS")
    if design["donors_per_superpopulation"] != 5:
        raise Stage1PreparationError("Stage 1 requires exactly five donors per superpopulation")
    if tuple(design["superpopulations"]) != REQUIRED_GROUPS:
        raise Stage1PreparationError("Stage 1 superpopulations must be AFR, AMR, EAS, EUR, SAS")
    if tuple(design["backgrounds"]) != ("defined", "complex"):
        raise Stage1PreparationError("Stage 1 requires defined and complex microbial backgrounds")
    if tuple(float(value) for value in design["spike_fractions"]) != REQUIRED_FRACTIONS:
        raise Stage1PreparationError("Stage 1 spike fractions must be 0.1%, 1%, 5%, and 10%")
    if set(design["analyses"]) != {"leave_one_donor_out", "leave_one_superpopulation_out"}:
        raise Stage1PreparationError("Stage 1 requires both leave-one-out analyses")
    sources = design["allowed_assembly_sources"]
    if not isinstance(sources, list) or not {"HPRC_R2", "HGSVC3"}.issubset(set(sources)):
        raise Stage1PreparationError("Stage 1 must permit HPRC_R2 and HGSVC3 augmentation")
    return int(design["donors_per_superpopulation"])


def prepare_stage1(
    *,
    design_path: Path,
    donors_path: Path,
    excluded_donors_path: Path,
    output: Path,
) -> dict[str, Any]:
    """Validate a complete, independent 25-donor Stage 1 manifest.

    The donor file is a public metadata manifest.  It must already contain a
    checksum-pinned assembly accession for each selected donor; this function
    intentionally refuses to infer or download one.
    """

    design = _load_design(design_path)
    per_group = _validate_design(design)
    donors = read_tsv(donors_path)
    excluded = read_tsv(excluded_donors_path)
    if not donors:
        raise Stage1PreparationError("Stage 1 donor manifest is empty")
    if set(DONOR_FIELDS).difference(donors[0]):
        raise Stage1PreparationError("Stage 1 donor manifest has invalid columns")
    if "donor_id" not in excluded[0]:
        raise Stage1PreparationError("excluded donor manifest has no donor_id column")
    excluded_ids = {row["donor_id"] for row in excluded}
    allowed_sources = set(design["allowed_assembly_sources"])
    selected: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in donors:
        donor = row["donor_id"]
        if donor in seen:
            raise Stage1PreparationError(f"duplicate Stage 1 donor {donor}")
        seen.add(donor)
        if donor in excluded_ids:
            raise Stage1PreparationError(f"Stage 1 donor overlaps Gate A panel: {donor}")
        if row["superpopulation"] not in REQUIRED_GROUPS:
            raise Stage1PreparationError(f"{donor}: invalid superpopulation")
        if row["assembly_source"] not in allowed_sources:
            raise Stage1PreparationError(f"{donor}: unsupported assembly source")
        if not row["assembly_accession"]:
            raise Stage1PreparationError(f"{donor}: missing assembly accession")
        try:
            if int(row["selection_rank"]) < 1:
                raise ValueError
        except ValueError as error:
            raise Stage1PreparationError(f"{donor}: invalid selection rank") from error
        selected.append({field: row[field] for field in DONOR_FIELDS})
    counts = Counter(row["superpopulation"] for row in selected)
    if any(counts[group] != per_group for group in REQUIRED_GROUPS):
        raise Stage1PreparationError(
            "Stage 1 requires five independent donors in every superpopulation; "
            f"observed={dict(sorted(counts.items()))}"
        )
    selected.sort(key=lambda row: (row["superpopulation"], int(row["selection_rank"]), row["donor_id"]))
    report = {
        "schema_version": 1,
        "status": "prepared_no_outcomes",
        "gate_a_status": design["gate_a_status"],
        "design_sha256": _sha256(design_path),
        "donor_manifest_sha256": sha256_bytes(canonical_tsv(selected, DONOR_FIELDS)),
        "excluded_panel_sha256": _sha256(excluded_donors_path),
        "donor_count": len(selected),
        "donors_per_superpopulation": per_group,
        "counts_by_superpopulation": dict(sorted(counts.items())),
        "assembly_sources": dict(sorted(Counter(row["assembly_source"] for row in selected).items())),
        "backgrounds": design["backgrounds"],
        "spike_fractions": design["spike_fractions"],
        "analyses": sorted(design["analyses"]),
        "outcomes_generated": False,
    }
    write_json_atomic(report, output)
    return report
