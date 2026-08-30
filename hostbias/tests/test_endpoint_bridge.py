import json
from pathlib import Path

import jsonschema
from typer.testing import CliRunner

from hostbias.cli import app
from hostbias.endpoint_bridge import aggregate_sample_endpoint


ROOT = Path(__file__).parents[1]
FIXTURES = Path(__file__).parent / "fixtures"


def test_endpoint_aggregate_is_privacy_safe_and_schema_valid() -> None:
    payload = aggregate_sample_endpoint(
        FIXTURES / "alignments.tsv",
        FIXTURES / "bins.tsv",
        FIXTURES / "bin_qc.tsv",
        "T01",
        "source",
        ROOT / "config" / "thresholds.yaml",
    )
    schema = json.loads(
        (ROOT / "schemas" / "sample_endpoint_aggregate.schema.json").read_text(
            encoding="utf-8"
        )
    )
    jsonschema.validate(payload, schema)
    assert payload["p_count"] == 1.0
    assert payload["endpoint_bin_human_tiers"]["dominant"] == 1
    serialized = json.dumps(payload)
    assert "human_clear" not in serialized
    assert "bin_novel" not in serialized
    assert str(FIXTURES) not in serialized


def test_bridge_cli_commands_emit_real_artifacts(tmp_path: Path) -> None:
    qc = tmp_path / "qc.json"
    result = CliRunner().invoke(
        app,
        [
            "assembly-qc",
            "--assembly",
            str(FIXTURES / "assembly.fa"),
            "--sample-id",
            "T01",
            "--filter-mode",
            "source",
            "--output",
            str(qc),
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(qc.read_text())["contig_count"] == 3

    endpoint = tmp_path / "endpoint.json"
    result = CliRunner().invoke(
        app,
        [
            "aggregate-endpoint",
            "--alignments",
            str(FIXTURES / "alignments.tsv"),
            "--contig-bins",
            str(FIXTURES / "bins.tsv"),
            "--bin-qc",
            str(FIXTURES / "bin_qc.tsv"),
            "--sample-id",
            "T01",
            "--filter-mode",
            "source",
            "--thresholds",
            str(ROOT / "config" / "thresholds.yaml"),
            "--output",
            str(endpoint),
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(endpoint.read_text())["p_bp"] == 1.0


def test_no_selected_bins_is_valid_zero_propagation(tmp_path: Path) -> None:
    bins = tmp_path / "bins.tsv"
    qc = tmp_path / "qc.tsv"
    bins.write_text("sample_id\tcontig_id\tbin_id\n", encoding="utf-8")
    qc.write_text(
        "sample_id\tbin_id\tdas_tool_selected\tcheckm2_completeness\t"
        "checkm2_contamination\tgunc_pass\tgtdb_domain\tgtdb_genus\t"
        "gtdb_species\n",
        encoding="utf-8",
    )
    payload = aggregate_sample_endpoint(
        FIXTURES / "alignments.tsv",
        bins,
        qc,
        "T01",
        "source",
    )
    assert payload["endpoint_bin_count"] == 0
    assert payload["propagated_human_contig_count"] == 0
    assert payload["p_count"] == 0
