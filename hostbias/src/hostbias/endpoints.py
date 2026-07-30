"""Quality-passing novel microbial bin and propagation endpoints."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

from hostbias.labeling import ContigCall, ContigLabel
from hostbias.schemas import BinQcRow, ContigBinRow, SchemaError, assert_unique


@dataclass(frozen=True)
class EndpointThresholds:
    min_completeness: float = 0.50
    max_contamination: float = 0.05
    material_human_fraction: float = 0.10
    dominant_human_fraction: float = 0.50


def is_endpoint_bin(row: BinQcRow, thresholds: EndpointThresholds | None = None) -> bool:
    """Apply the preregistered endpoint definition exactly."""

    thresholds = thresholds or EndpointThresholds()
    named_genus = bool(row.gtdb_genus and row.gtdb_genus.strip() not in {"g__", "NA"})
    no_species = not row.gtdb_species or row.gtdb_species.strip() in {"s__", "NA"}
    return (
        row.das_tool_selected
        and row.checkm2_completeness >= thresholds.min_completeness
        and row.checkm2_contamination < thresholds.max_contamination
        and row.gunc_pass
        and row.gtdb_domain in {"Bacteria", "Archaea"}
        and named_genus
        and no_species
    )


@dataclass(frozen=True)
class SampleEndpoint:
    sample_id: str
    human_contig_count: int
    propagated_human_contig_count: int
    human_bp: int
    propagated_human_bp: int
    p_count: float
    p_bp: float
    denominator_state: str
    endpoint_bin_count: int
    endpoint_bins_with_human: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class BinHumanFraction:
    sample_id: str
    bin_id: str
    total_bp: int
    human_bp: int
    human_fraction: float
    tier: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def calculate_endpoints(
    calls: Iterable[ContigCall],
    contig_bins: Iterable[ContigBinRow],
    bin_qc: Iterable[BinQcRow],
    thresholds: EndpointThresholds | None = None,
    expected_samples: Iterable[str] | None = None,
) -> tuple[list[SampleEndpoint], list[BinHumanFraction]]:
    """Join calls to bins without permitting duplicate-key inflation."""

    thresholds = thresholds or EndpointThresholds()
    calls = list(calls)
    contig_bins = list(contig_bins)
    bin_qc = list(bin_qc)
    assert_unique(calls, ("sample_id", "contig_id"))
    assert_unique(contig_bins, ("sample_id", "contig_id"))
    assert_unique(bin_qc, ("sample_id", "bin_id"))

    call_by_key = {(row.sample_id, row.contig_id): row for row in calls}
    bin_by_contig = {
        (row.sample_id, row.contig_id): row.bin_id for row in contig_bins
    }
    endpoint_bins = {
        (row.sample_id, row.bin_id)
        for row in bin_qc
        if is_endpoint_bin(row, thresholds)
    }
    sample_ids = (
        set(expected_samples)
        if expected_samples is not None
        else {row.sample_id for row in calls} | {row.sample_id for row in bin_qc}
    )

    bin_lengths: dict[tuple[str, str], int] = {}
    bin_human_lengths: dict[tuple[str, str], int] = {}
    for contig_key, bin_id in bin_by_contig.items():
        call = call_by_key.get(contig_key)
        if call is None:
            raise SchemaError(f"bin mapping has no contig call: {contig_key}")
        bin_key = (call.sample_id, bin_id)
        bin_lengths[bin_key] = bin_lengths.get(bin_key, 0) + call.contig_length
        if call.label == ContigLabel.HUMAN:
            bin_human_lengths[bin_key] = (
                bin_human_lengths.get(bin_key, 0) + call.contig_length
            )

    fraction_rows: list[BinHumanFraction] = []
    for bin_key in sorted(endpoint_bins):
        total_bp = bin_lengths.get(bin_key, 0)
        human_bp = bin_human_lengths.get(bin_key, 0)
        fraction = human_bp / total_bp if total_bp else 0.0
        if fraction >= thresholds.dominant_human_fraction:
            tier = "dominant"
        elif fraction >= thresholds.material_human_fraction:
            tier = "material"
        elif human_bp:
            tier = "contact"
        else:
            tier = "none"
        fraction_rows.append(
            BinHumanFraction(*bin_key, total_bp, human_bp, fraction, tier)
        )

    sample_rows: list[SampleEndpoint] = []
    for sample_id in sorted(sample_ids):
        human_calls = [
            row
            for row in calls
            if row.sample_id == sample_id and row.label == ContigLabel.HUMAN
        ]
        propagated = [
            row
            for row in human_calls
            if (sample_id, bin_by_contig.get((sample_id, row.contig_id), ""))
            in endpoint_bins
        ]
        human_bp = sum(row.contig_length for row in human_calls)
        propagated_bp = sum(row.contig_length for row in propagated)
        endpoint_for_sample = {
            bin_id for current_sample, bin_id in endpoint_bins if current_sample == sample_id
        }
        with_human = {
            bin_by_contig[(sample_id, row.contig_id)]
            for row in propagated
        }
        sample_rows.append(
            SampleEndpoint(
                sample_id=sample_id,
                human_contig_count=len(human_calls),
                propagated_human_contig_count=len(propagated),
                human_bp=human_bp,
                propagated_human_bp=propagated_bp,
                p_count=len(propagated) / len(human_calls) if human_calls else 0.0,
                p_bp=propagated_bp / human_bp if human_bp else 0.0,
                denominator_state="observed" if human_calls else "zero",
                endpoint_bin_count=len(endpoint_for_sample),
                endpoint_bins_with_human=len(with_human),
            )
        )
    return sample_rows, fraction_rows
