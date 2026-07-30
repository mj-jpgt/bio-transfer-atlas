"""Configuration and manifest validation."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import jsonschema
import yaml


class ValidationError(ValueError):
    """Raised when a Hostbias input contract is invalid."""


@dataclass(frozen=True)
class ValidatedInputs:
    """Validated configuration plus resolved sample rows."""

    config_path: Path
    root: Path
    config: dict[str, Any]
    manifest_path: Path
    samples: tuple[dict[str, str], ...]


def _load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _validate(instance: Any, schema: dict[str, Any], label: str) -> None:
    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(instance), key=lambda error: list(error.path))
    if errors:
        details = "; ".join(
            f"{'.'.join(map(str, error.path)) or '<root>'}: {error.message}" for error in errors
        )
        raise ValidationError(f"{label} failed schema validation: {details}")


def _resolve_root(config_path: Path) -> Path:
    """Resolve the project root from a config path under ``config/``."""

    config_path = config_path.resolve()
    if config_path.parent.name == "config":
        return config_path.parent.parent
    return config_path.parent


def load_and_validate(config_path: str | Path) -> ValidatedInputs:
    config_path = Path(config_path).resolve()
    if not config_path.is_file():
        raise ValidationError(f"configuration does not exist: {config_path}")
    root = _resolve_root(config_path)
    schema_dir = root / "schemas"

    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValidationError("configuration must be a YAML mapping")
    _validate(config, _load_json(schema_dir / "config.schema.json"), "configuration")

    manifest_path = Path(config["paths"]["sample_manifest"])
    if not manifest_path.is_absolute():
        manifest_path = root / manifest_path
    if not manifest_path.is_file():
        raise ValidationError(f"sample manifest does not exist: {manifest_path}")

    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = tuple(dict(row) for row in reader)
        fieldnames = set(reader.fieldnames or ())

    row_schema = _load_json(schema_dir / "sample_manifest.schema.json")
    expected = set(row_schema["required"])
    if fieldnames != expected:
        missing = sorted(expected - fieldnames)
        extra = sorted(fieldnames - expected)
        raise ValidationError(f"sample manifest columns differ; missing={missing}, extra={extra}")
    if not rows:
        raise ValidationError("sample manifest contains no samples")
    for index, row in enumerate(rows, start=2):
        _validate(row, row_schema, f"sample manifest line {index}")

    sample_ids = [row["sample_id"] for row in rows]
    accessions = [row["accession"] for row in rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValidationError("sample_id values must be unique")
    if len(accessions) != len(set(accessions)):
        raise ValidationError("accession values must be unique")
    for row in rows:
        if row["fastq1_url"] == row["fastq2_url"]:
            raise ValidationError(f"{row['sample_id']}: mate URLs must differ")
        if row["fastq1_md5"].lower() == row["fastq2_md5"].lower():
            # Zero placeholders are allowed only in the distributed example.
            placeholder = "0" * 32
            if row["fastq1_md5"] != placeholder or manifest_path.name.endswith(".example.tsv") is False:
                raise ValidationError(f"{row['sample_id']}: mate checksums must differ")

    return ValidatedInputs(
        config_path=config_path,
        root=root,
        config=config,
        manifest_path=manifest_path.resolve(),
        samples=rows,
    )
