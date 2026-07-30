import json
from pathlib import Path

from typer.testing import CliRunner

from hostbias.cli import app, run_analysis


FIXTURES = Path(__file__).parent / "fixtures"


def test_end_to_end_analysis_writes_all_audit_artifacts(tmp_path: Path) -> None:
    status = run_analysis(
        FIXTURES / "alignments.tsv",
        FIXTURES / "bins.tsv",
        FIXTURES / "bin_qc.tsv",
        FIXTURES / "sample_groups.tsv",
        FIXTURES / "control_alignments.tsv",
        FIXTURES / "control_truth.tsv",
        FIXTURES / "sensitivities.tsv",
        tmp_path,
        FIXTURES / "thresholds.yaml",
    )
    assert status == "FAIL"
    expected = {
        "effective_thresholds.json",
        "contig_calls.json",
        "sample_endpoints.json",
        "bin_human_fractions.json",
        "control_results.json",
        "analysis_statistics.json",
        "GATE_A_VERDICT.json",
        "GATE_A_VERDICT.md",
    }
    assert expected == {path.name for path in tmp_path.iterdir()}
    payload = json.loads((tmp_path / "GATE_A_VERDICT.json").read_text())
    assert payload["status"] == "FAIL"
    assert payload["first_failed_criterion"] == "ratio_ci_excludes_one"


def test_cli_help_exposes_single_complete_entrypoint() -> None:
    result = CliRunner().invoke(app, ["analyze", "--help"])
    assert result.exit_code == 0
    assert "--control-alignments" in result.stdout
    assert "--sensitivities" in result.stdout
