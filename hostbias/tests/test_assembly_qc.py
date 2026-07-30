import json
from pathlib import Path

import jsonschema
import pytest

from hostbias.assembly_qc import assembly_qc, iter_fasta
from hostbias.schemas import SchemaError


ROOT = Path(__file__).parents[1]
FIXTURES = Path(__file__).parent / "fixtures"


def test_assembly_qc_is_aggregate_and_schema_valid() -> None:
    result = assembly_qc(FIXTURES / "assembly.fa", "T01", "source")
    schema = json.loads(
        (ROOT / "schemas" / "assembly_qc.schema.json").read_text(encoding="utf-8")
    )
    jsonschema.validate(result, schema)
    assert result["contig_count"] == 3
    assert result["total_bp"] == 22
    assert result["n50_bp"] == 12
    assert result["n90_bp"] == 4
    serialized = json.dumps(result)
    assert "c1" not in serialized
    assert str(FIXTURES) not in serialized


def test_fasta_parser_rejects_duplicate_identifiers(tmp_path: Path) -> None:
    path = tmp_path / "duplicate.fa"
    path.write_text(">same\nAAAA\n>same\nCCCC\n", encoding="utf-8")
    with pytest.raises(SchemaError, match="duplicate FASTA identifier"):
        list(iter_fasta(path))


def test_fasta_parser_rejects_non_dna_symbols(tmp_path: Path) -> None:
    path = tmp_path / "invalid.fa"
    path.write_text(">c\nACGTZ\n", encoding="utf-8")
    with pytest.raises(SchemaError, match="invalid DNA symbols"):
        assembly_qc(path, "T01", "source")
