from __future__ import annotations

import json
from pathlib import Path

import pytest

from hostbias.stage1_prep import Stage1PreparationError, prepare_stage1


def _write_design(path: Path) -> None:
    path.write_text(
        """gate_a_status: PASS_REQUIRED
donors_per_superpopulation: 5
superpopulations: [AFR, AMR, EAS, EUR, SAS]
backgrounds: [defined, complex]
spike_fractions: [0.001, 0.01, 0.05, 0.1]
analyses: [leave_one_donor_out, leave_one_superpopulation_out]
allowed_assembly_sources: [HPRC_R2, HGSVC3]
""",
        encoding="utf-8",
    )


def _write_donors(path: Path) -> None:
    rows = ["donor_id\tsuperpopulation\tpopulation_code\tassembly_source\tassembly_accession\tselection_rank"]
    for group in ("AFR", "AMR", "EAS", "EUR", "SAS"):
        for rank in range(1, 6):
            source = "HGSVC3" if group == "AMR" else "HPRC_R2"
            rows.append(f"{group}{rank}\t{group}\tPOP\t{source}\tGCA_{group}{rank}\t{rank}")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def test_stage1_preparation_writes_aggregate_only_checkpoint(tmp_path: Path) -> None:
    design, donors, excluded, output = (tmp_path / name for name in ("design.yaml", "donors.tsv", "excluded.tsv", "stage1.json"))
    _write_design(design)
    _write_donors(donors)
    excluded.write_text("donor_id\nGATEA01\n", encoding="utf-8")

    report = prepare_stage1(design_path=design, donors_path=donors, excluded_donors_path=excluded, output=output)

    assert report["status"] == "prepared_no_outcomes"
    assert report["donor_count"] == 25
    assert report["counts_by_superpopulation"] == {group: 5 for group in ("AFR", "AMR", "EAS", "EUR", "SAS")}
    assert json.loads(output.read_text())["outcomes_generated"] is False


def test_stage1_preparation_rejects_gate_a_panel_overlap(tmp_path: Path) -> None:
    design, donors, excluded = tmp_path / "design.yaml", tmp_path / "donors.tsv", tmp_path / "excluded.tsv"
    _write_design(design)
    _write_donors(donors)
    excluded.write_text("donor_id\nAFR1\n", encoding="utf-8")

    with pytest.raises(Stage1PreparationError, match="overlaps Gate A panel"):
        prepare_stage1(design_path=design, donors_path=donors, excluded_donors_path=excluded, output=tmp_path / "out.json")


def test_stage1_preparation_rejects_incomplete_group(tmp_path: Path) -> None:
    design, donors, excluded = tmp_path / "design.yaml", tmp_path / "donors.tsv", tmp_path / "excluded.tsv"
    _write_design(design)
    _write_donors(donors)
    rows = donors.read_text().splitlines()
    donors.write_text("\n".join(row for row in rows if not row.startswith("AMR5\t")) + "\n", encoding="utf-8")
    excluded.write_text("donor_id\nGATEA01\n", encoding="utf-8")

    with pytest.raises(Stage1PreparationError, match="five independent donors"):
        prepare_stage1(design_path=design, donors_path=donors, excluded_donors_path=excluded, output=tmp_path / "out.json")
