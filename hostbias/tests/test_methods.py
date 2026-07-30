from pathlib import Path

import pytest

from hostbias.controls import evaluate_controls
from hostbias.endpoints import calculate_endpoints, is_endpoint_bin
from hostbias.labeling import ContigCall, ContigLabel, LabelThresholds, label_contigs
from hostbias.schemas import (
    AlignmentRow,
    BinQcRow,
    ContigBinRow,
    ControlTruthRow,
    SchemaError,
    read_tsv,
)


FIXTURES = Path(__file__).parent / "fixtures"


def test_competitive_labeling_and_identity_sensitivity() -> None:
    rows = read_tsv(FIXTURES / "alignments.tsv", AlignmentRow)
    calls = {row.contig_id: row for row in label_contigs(rows)}
    assert calls["human_clear"].label == ContigLabel.HUMAN
    assert calls["ambiguous"].label == ContigLabel.AMBIGUOUS
    assert calls["microbial"].label == ContigLabel.NON_HUMAN
    calls_98 = {row.contig_id: row for row in label_contigs(
        rows, LabelThresholds(min_identity=0.98)
    )}
    assert calls_98["human_clear"].label == ContigLabel.HUMAN


def test_no_hit_placeholder_is_non_human() -> None:
    row = AlignmentRow("T", "unmatched", 1000, "none", 0, 0, 0, 0, 0)
    assert label_contigs([row])[0].label == ContigLabel.NON_HUMAN
    with pytest.raises(SchemaError, match="zero-valued"):
        AlignmentRow("T", "bad", 1000, "none", 1, 0, 0, 0, 0)


def test_endpoint_definition_and_sample_rates() -> None:
    calls = label_contigs(read_tsv(FIXTURES / "alignments.tsv", AlignmentRow))
    bins = read_tsv(FIXTURES / "bins.tsv", ContigBinRow)
    qc = read_tsv(FIXTURES / "bin_qc.tsv", BinQcRow)
    assert is_endpoint_bin(qc[0])
    assert not is_endpoint_bin(qc[1])
    samples, fractions = calculate_endpoints(
        calls, bins, qc, expected_samples=["T01", "T02"]
    )
    t01, t02 = samples
    assert (t01.p_count, t01.p_bp) == (1.0, 1.0)
    assert t01.endpoint_bins_with_human == 1
    assert fractions[0].human_fraction == pytest.approx(0.5)
    assert fractions[0].tier == "dominant"
    assert t02.denominator_state == "zero"
    assert t02.p_count == 0


def test_join_rejects_unlabelled_bin_contig() -> None:
    calls = [ContigCall("S", "c1", 100, ContigLabel.HUMAN, 100.0, None, None)]
    bins = [ContigBinRow("S", "missing", "b1")]
    qc = [BinQcRow("S", "b1", True, 0.8, 0.01, True, "Bacteria", "g__X", None)]
    with pytest.raises(SchemaError, match="has no contig call"):
        calculate_endpoints(calls, bins, qc)


def test_control_evaluation_uses_bp_specificity() -> None:
    calls = [
        ContigCall("C", "h1", 1000, ContigLabel.HUMAN, 10, None, None),
        ContigCall("C", "h2", 1000, ContigLabel.HUMAN, 10, None, None),
        ContigCall("C", "m1", 1, ContigLabel.HUMAN, 10, None, None),
        ContigCall("C", "m2", 9999, ContigLabel.NON_HUMAN, None, 10, None),
    ]
    truth = [
        ControlTruthRow("C", "h1", "human", 1000),
        ControlTruthRow("C", "h2", "human", 1000),
        ControlTruthRow("C", "m1", "microbial", 1),
        ControlTruthRow("C", "m2", "microbial", 9999),
    ]
    result = evaluate_controls(calls, truth)
    assert result.human_sensitivity == 1
    assert result.microbial_false_positive_rate == 0.5
    assert result.microbial_false_positive_bp_rate == 0.0001
    assert result.passed
