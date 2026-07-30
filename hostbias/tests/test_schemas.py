from pathlib import Path

import pytest

from hostbias.schemas import (
    AlignmentRow,
    BinQcRow,
    ContigBinRow,
    SchemaError,
    assert_unique,
    read_tsv,
)


FIXTURES = Path(__file__).parent / "fixtures"


def test_reads_strict_tables() -> None:
    alignments = read_tsv(FIXTURES / "alignments.tsv", AlignmentRow)
    bins = read_tsv(FIXTURES / "bins.tsv", ContigBinRow)
    qc = read_tsv(FIXTURES / "bin_qc.tsv", BinQcRow)
    assert len(alignments) == 5
    assert len(bins) == 3
    assert qc[0].gtdb_species is None


def test_rejects_wrong_header(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("sample_id\tcontig_id\nT01\tc1\n", encoding="utf-8")
    with pytest.raises(SchemaError, match="expected header"):
        read_tsv(path, AlignmentRow)


def test_rejects_duplicate_natural_key() -> None:
    row = ContigBinRow("T01", "c1", "b1")
    with pytest.raises(SchemaError, match="duplicate key"):
        assert_unique([row, row], ("sample_id", "contig_id"))
