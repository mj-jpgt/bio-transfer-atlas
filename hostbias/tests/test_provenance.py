from __future__ import annotations

import json
from pathlib import Path

from hostbias.config import load_and_validate
from hostbias.provenance import build_provenance, write_json_atomic


PROJECT = Path(__file__).parents[1]


def test_provenance_is_aggregate_safe(tmp_path: Path) -> None:
    inputs = load_and_validate(PROJECT / "config" / "config.example.yaml")
    provenance = build_provenance(inputs)
    serialized = json.dumps(provenance)

    assert provenance["inputs"]["sample_manifest"]["sample_count"] == 2
    assert provenance["privacy"]["environment_variables_recorded"] is False
    assert str(PROJECT.resolve()) not in serialized
    assert "fastq1_url" not in serialized

    output = tmp_path / "provenance.json"
    write_json_atomic(provenance, output)
    assert json.loads(output.read_text()) == provenance
