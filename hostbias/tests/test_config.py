from __future__ import annotations

import csv
from pathlib import Path

import pytest
import yaml

from hostbias.config import ValidationError, load_and_validate


PROJECT = Path(__file__).parents[1]


def test_example_configuration_is_valid() -> None:
    inputs = load_and_validate(PROJECT / "config" / "config.example.yaml")
    assert inputs.config["experiment"]["seed"] == 20260729
    assert {sample["cohort"] for sample in inputs.samples} == {"Tanzania", "Netherlands"}
    assert inputs.config["databases"]["gtdb_release"] == "R220"
    assert inputs.config["binning"]["maxbin_markerset"] == 40


def test_duplicate_sample_ids_are_rejected(tmp_path: Path) -> None:
    config = yaml.safe_load((PROJECT / "config" / "config.example.yaml").read_text())
    source = PROJECT / "config" / "stage0_samples.example.tsv"
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows[1]["sample_id"] = rows[0]["sample_id"]
    manifest = tmp_path / "samples.tsv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    config["paths"]["sample_manifest"] = str(manifest)
    config_path = PROJECT / "config" / "_duplicate.test.yaml"
    try:
        config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
        with pytest.raises(ValidationError, match="sample_id values must be unique"):
            load_and_validate(config_path)
    finally:
        config_path.unlink(missing_ok=True)


def test_unknown_configuration_key_is_rejected(tmp_path: Path) -> None:
    config = yaml.safe_load((PROJECT / "config" / "config.example.yaml").read_text())
    config["surprise"] = True
    config_path = PROJECT / "config" / "_invalid.test.yaml"
    try:
        config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
        with pytest.raises(ValidationError, match="Additional properties"):
            load_and_validate(config_path)
    finally:
        config_path.unlink(missing_ok=True)
