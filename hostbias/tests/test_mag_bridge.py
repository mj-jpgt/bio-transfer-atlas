from pathlib import Path

import pytest
from typer.testing import CliRunner

from hostbias.cli import app
from hostbias.mag_bridge import bins_to_scaffolds2bin, depth_to_maxbin_abundance
from hostbias.schemas import SchemaError


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
