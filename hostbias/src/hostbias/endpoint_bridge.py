"""Privacy-safe endpoint aggregation for one assembly/filter-mode unit."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

from hostbias.endpoints import EndpointThresholds, calculate_endpoints
from hostbias.labeling import LabelThresholds, label_contigs
from hostbias.provenance import sha256_file
from hostbias.schemas import (
    AlignmentRow,
    BinQcRow,
    ContigBinRow,
    SchemaError,
    read_tsv,
)


def load_preregistered_method_thresholds(
    path: Path | None,
) -> tuple[LabelThresholds, EndpointThresholds]:
    """Load the frozen preregistration shape used by ``config/thresholds.yaml``."""

    if path is None:
        return LabelThresholds(), EndpointThresholds()
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise SchemaError("thresholds must be a YAML mapping")
    labels = raw.get("human_contig_label")
    endpoint = raw.get("endpoint_bin")
    if not isinstance(labels, dict) or not isinstance(endpoint, dict):
        raise SchemaError(
            "thresholds require human_contig_label and endpoint_bin mappings"
        )
    tiers = endpoint.get("human_fraction_tiers")
    if not isinstance(tiers, dict):
        raise SchemaError("endpoint_bin.human_fraction_tiers must be a mapping")
    try:
        return (
            LabelThresholds(
                min_aligned_bp=int(labels["minimum_aligned_bp"]),
                min_identity=float(labels["primary_identity_fraction"]),
                min_query_coverage=float(
                    labels["minimum_query_coverage_fraction"]
                ),
                min_mapq=float(labels["minimum_mapping_quality"]),
                min_human_score_ratio=float(
                    labels["minimum_human_to_gtdb_score_ratio"]
                ),
            ),
            EndpointThresholds(
                min_completeness=float(
                    endpoint["minimum_checkm2_completeness_fraction"]
                ),
                max_contamination=float(
                    endpoint["maximum_checkm2_contamination_fraction"]
                ),
                material_human_fraction=float(tiers["material"]),
                dominant_human_fraction=float(tiers["dominant"]),
            ),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise SchemaError(f"invalid preregistered method thresholds: {exc}") from exc


def aggregate_sample_endpoint(
    alignments_path: Path,
    contig_bins_path: Path,
    bin_qc_path: Path,
    sample_id: str,
    filter_mode: str,
    thresholds_path: Path | None = None,
) -> dict[str, Any]:
    """Calculate one aggregate endpoint record and omit all contig/bin identifiers."""

    if filter_mode not in {"source", "strict"}:
        raise SchemaError("filter_mode must be 'source' or 'strict'")
    alignments = read_tsv(alignments_path, AlignmentRow)
    contig_bins = read_tsv(contig_bins_path, ContigBinRow)
    bin_qc = read_tsv(bin_qc_path, BinQcRow)
    observed_samples = (
        {row.sample_id for row in alignments}
        | {row.sample_id for row in contig_bins}
        | {row.sample_id for row in bin_qc}
    )
    if observed_samples != {sample_id}:
        raise SchemaError(
            f"endpoint inputs must contain only sample {sample_id!r}; "
            f"found {sorted(observed_samples)}"
        )
    label_thresholds, endpoint_thresholds = load_preregistered_method_thresholds(
        thresholds_path
    )
    calls = label_contigs(alignments, label_thresholds)
    endpoints, fractions = calculate_endpoints(
        calls,
        contig_bins,
        bin_qc,
        endpoint_thresholds,
        expected_samples=[sample_id],
    )
    endpoint = endpoints[0]
    tiers = Counter(row.tier for row in fractions)
    return {
        "schema_version": "1.0",
        "sample_id": sample_id,
        "filter_mode": filter_mode,
        "input_sha256": {
            "alignments": sha256_file(alignments_path),
            "contig_bins": sha256_file(contig_bins_path),
            "bin_qc": sha256_file(bin_qc_path),
        },
        "human_contig_count": endpoint.human_contig_count,
        "propagated_human_contig_count": endpoint.propagated_human_contig_count,
        "human_bp": endpoint.human_bp,
        "propagated_human_bp": endpoint.propagated_human_bp,
        "p_count": endpoint.p_count,
        "p_bp": endpoint.p_bp,
        "denominator_state": endpoint.denominator_state,
        "endpoint_bin_count": endpoint.endpoint_bin_count,
        "endpoint_bins_with_human": endpoint.endpoint_bins_with_human,
        "endpoint_bin_human_tiers": {
            tier: tiers.get(tier, 0)
            for tier in ("none", "contact", "material", "dominant")
        },
        "maximum_endpoint_bin_human_fraction": max(
            (row.human_fraction for row in fractions), default=0.0
        ),
        "privacy": {
            "contains_sequences": False,
            "contains_contig_identifiers": False,
            "contains_bin_identifiers": False,
            "contains_filesystem_paths": False,
        },
    }
