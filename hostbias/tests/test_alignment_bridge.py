import json
from pathlib import Path

import jsonschema
import pytest

from hostbias.alignment_bridge import (
    build_alignment_table,
    load_alignment_run_spec,
    run_alignment_bridge,
)
from hostbias.schemas import AlignmentRow, SchemaError, read_tsv


ROOT = Path(__file__).parents[1]
FIXTURES = Path(__file__).parent / "fixtures"


def _spec(tmp_path: Path) -> Path:
    assembly = FIXTURES / "assembly.fa"
    import hashlib

    payload = {
        "schema_version": "1.0",
        "sample_id": "T01",
        "filter_mode": "source",
        "assembly": {
            "path": str(assembly),
            "sha256": hashlib.sha256(assembly.read_bytes()).hexdigest(),
        },
        "mapper": {"version": "2.28", "preset": "asm20"},
        "references": [
            {
                "domain": "human",
                "reference_id": "hprc-balanced-v1",
                "reference_sha256": "a" * 64,
                "paf_path": str(FIXTURES / "human.paf"),
            },
            {
                "domain": "gtdb",
                "reference_id": "gtdb-r220",
                "reference_sha256": "b" * 64,
                "paf_path": str(FIXTURES / "gtdb.paf"),
            },
        ],
    }
    path = tmp_path / "spec.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_paf_bridge_builds_exact_analysis_contract(tmp_path: Path) -> None:
    spec_path = _spec(tmp_path)
    output_tsv = tmp_path / "alignments.tsv"
    output_manifest = tmp_path / "manifest.json"
    run_alignment_bridge(spec_path, output_tsv, output_manifest)

    rows = read_tsv(output_tsv, AlignmentRow)
    assert len(rows) == 4
    assert {(row.contig_id, row.target_domain) for row in rows} == {
        ("c1", "human"),
        ("c1", "gtdb"),
        ("c2", "human"),
        ("c3", "none"),
    }
    manifest = json.loads(output_manifest.read_text(encoding="utf-8"))
    schema = json.loads(
        (ROOT / "schemas" / "competitive_alignment_manifest.schema.json").read_text(
            encoding="utf-8"
        )
    )
    jsonschema.validate(manifest, schema)
    serialized = json.dumps(manifest)
    assert "c1" not in serialized
    assert str(FIXTURES) not in serialized
    assert manifest["aggregate"]["no_hit_contig_count"] == 1


def test_paf_query_length_must_match_assembly(tmp_path: Path) -> None:
    spec = load_alignment_run_spec(_spec(tmp_path))
    bad_paf = tmp_path / "bad.paf"
    bad_paf.write_text(
        "c1\t13\t0\t12\t+\tt\t100\t0\t12\t12\t12\t60\tAS:i:1\n",
        encoding="utf-8",
    )
    bad_source = type(spec.sources[0])(
        "human", "bad", "c" * 64, bad_paf
    )
    bad_spec = type(spec)(
        spec.sample_id,
        spec.filter_mode,
        spec.assembly_path,
        spec.assembly_sha256,
        spec.mapper_version,
        spec.mapper_preset,
        (bad_source, spec.sources[1]),
    )
    with pytest.raises(SchemaError, match="query length mismatch"):
        build_alignment_table(bad_spec)
