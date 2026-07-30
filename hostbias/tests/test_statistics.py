import json

from hostbias.controls import ControlResult
from hostbias.schemas import SensitivityRow
from hostbias.statistics import GroupedEndpoint, analyze
from hostbias.verdict import make_verdict, write_verdict


def _passing_rows() -> list[GroupedEndpoint]:
    rows = []
    for index in range(20):
        rows.append(
            GroupedEndpoint(
                f"T{index:02}", "tanzania", 0.10 + index / 1000, 0.08, 2
            )
        )
        rows.append(
            GroupedEndpoint(
                f"N{index:02}", "netherlands", 0.01 + index / 10000, 0.005, 1
            )
        )
    return rows


def _sensitivities(direction: bool = True) -> list[SensitivityRow]:
    return [
        SensitivityRow(
            analysis_id,
            "p_count",
            0.08 if direction else 0.01,
            0.01 if direction else 0.08,
        )
        for analysis_id in (
            "strict_pair",
            "identity_0.90",
            "identity_0.95",
            "identity_0.98",
        )
    ]


def _passing_controls() -> ControlResult:
    return ControlResult(100, 100, 0.98, 0.98, 0.0, 0.0, True)


def test_statistics_are_deterministic_and_gate_passes(tmp_path) -> None:
    first = analyze(_passing_rows(), 500, 1_000)
    second = analyze(_passing_rows(), 500, 1_000)
    assert first == second
    verdict = make_verdict(first, _passing_controls(), _sensitivities())
    assert verdict.status == "PASS"
    json_path, markdown_path = write_verdict(verdict, tmp_path)
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["pass_gate_a"] is True
    assert "Status:** PASS" in markdown_path.read_text(encoding="utf-8")


def test_scientific_failure_names_first_failed_link() -> None:
    rows = [
        GroupedEndpoint(
            row.sample_id,
            row.cohort,
            0.0,
            row.p_bp,
            0,
        )
        for row in _passing_rows()
    ]
    stats = analyze(rows, 200, 200)
    verdict = make_verdict(stats, _passing_controls(), _sensitivities())
    assert verdict.status == "FAIL"
    assert verdict.pass_gate_a is False
    assert verdict.first_failed_criterion == "minimum_tanzania_propagation"


def test_control_failure_is_operational_not_scientific() -> None:
    controls = ControlResult(10, 10, 0.5, 0.5, 0.0, 0.0, False)
    verdict = make_verdict(
        analyze(_passing_rows(), 200, 200), controls, _sensitivities()
    )
    assert verdict.status == "OPERATIONAL_FAILURE"
    assert verdict.pass_gate_a is None
    assert verdict.first_failed_criterion == "controls_pass"


def test_missing_sensitivity_is_operational_failure() -> None:
    verdict = make_verdict(
        analyze(_passing_rows(), 200, 200),
        _passing_controls(),
        _sensitivities()[:-1],
    )
    assert verdict.status == "OPERATIONAL_FAILURE"
    assert verdict.first_failed_criterion == "complete_sensitivity_matrix"
