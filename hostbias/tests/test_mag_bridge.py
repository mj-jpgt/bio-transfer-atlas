from pathlib import Path

import pytest
from typer.testing import CliRunner

from hostbias.cli import app
from hostbias.mag_bridge import (
    bins_to_scaffolds2bin,
    build_mag_contracts,
    depth_to_maxbin_abundance,
)
from hostbias.schemas import BinQcRow, ContigBinRow, SchemaError, read_tsv


FIXTURES = Path(__file__).parent / "fixtures"


def test_depth_translation_is_deterministic_and_headerless(tmp_path: Path) -> None:
    output = tmp_path / "abundance.tsv"
    assert depth_to_maxbin_abundance(FIXTURES / "contig_depth.tsv", output) == 3
    assert output.read_text(encoding="utf-8").splitlines() == [
        "c1\t10.25",
        "c2\t0",
        "c3\t3.5",
    ]


def test_bin_directory_translation_is_deterministic(tmp_path: Path) -> None:
    output = tmp_path / "scaffolds2bin.tsv"
    count = bins_to_scaffolds2bin(
        FIXTURES / "mag_bins", output, "metabat."
    )
    assert count == 3
    assert output.read_text(encoding="utf-8").splitlines() == [
        "c1\tmetabat.bin.1",
        "c2\tmetabat.bin.1",
        "c3\tmetabat.bin.2",
    ]


def test_duplicate_contig_across_bins_is_rejected(tmp_path: Path) -> None:
    bins = tmp_path / "bins"
    bins.mkdir()
    (bins / "one.fa").write_text(">same\nAAAA\n", encoding="utf-8")
    (bins / "two.fa").write_text(">same\nCCCC\n", encoding="utf-8")
    with pytest.raises(SchemaError, match="multiple bins"):
        bins_to_scaffolds2bin(bins, tmp_path / "out.tsv", "x.")


def test_mag_translation_commands_are_exposed(tmp_path: Path) -> None:
    output = tmp_path / "abundance.tsv"
    result = CliRunner().invoke(
        app,
        [
            "maxbin-abundance",
            "--depth",
            str(FIXTURES / "contig_depth.tsv"),
            "--output",
            str(output),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "3 contigs" in result.output


def test_mag_reports_join_into_exact_private_contracts(tmp_path: Path) -> None:
    contig_bins = tmp_path / "contig_bins.tsv"
    bin_qc = tmp_path / "bin_qc.tsv"
    count = build_mag_contracts(
        "T01",
        FIXTURES / "dastool_scaffolds2bin.tsv",
        FIXTURES / "checkm2_quality.tsv",
        FIXTURES / "gunc_maxcss.tsv",
        (
            FIXTURES / "gtdb_bac120.summary.tsv",
            FIXTURES / "gtdb_ar53.summary.tsv",
        ),
        contig_bins,
        bin_qc,
    )
    assert count == 2
    mappings = read_tsv(contig_bins, ContigBinRow)
    qc = read_tsv(bin_qc, BinQcRow)
    assert len(mappings) == 3
    assert qc[0].checkm2_completeness == 0.8
    assert qc[0].checkm2_contamination == 0.02
    assert qc[0].gtdb_genus == "g__Novel"
    assert qc[0].gtdb_species is None
    assert qc[1].gunc_pass is False


def test_mag_contract_rejects_missing_tool_bin(tmp_path: Path) -> None:
    gunc = tmp_path / "gunc.tsv"
    gunc.write_text(
        "genome\tpass.GUNC\nbin.1\tTrue\n",
        encoding="utf-8",
    )
    with pytest.raises(SchemaError, match="GUNC bins differ"):
        build_mag_contracts(
            "T01",
            FIXTURES / "dastool_scaffolds2bin.tsv",
            FIXTURES / "checkm2_quality.tsv",
            gunc,
            (
                FIXTURES / "gtdb_bac120.summary.tsv",
                FIXTURES / "gtdb_ar53.summary.tsv",
            ),
            tmp_path / "bins.tsv",
            tmp_path / "qc.tsv",
        )
