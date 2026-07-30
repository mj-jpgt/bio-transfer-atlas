"""Strict threshold configuration loading."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any, TypeVar

import jsonschema
import yaml  # type: ignore[import-untyped]

from hostbias.endpoints import EndpointThresholds
from hostbias.labeling import LabelThresholds
from hostbias.schemas import SchemaError
from hostbias.verdict import GateThresholds


T = TypeVar("T")


class ValidationError(ValueError):
    """Raised when a Hostbias workflow input contract is invalid."""


@dataclass(frozen=True)
class ValidatedInputs:
    """Validated workflow configuration plus resolved sample rows."""

    config_path: Path
    root: Path
    config: dict[str, Any]
    manifest_path: Path
    samples: tuple[dict[str, str], ...]


def _construct(section: str, cls: type[T], values: object) -> T:
    if values is None:
        return cls()
    if not isinstance(values, dict):
        raise SchemaError(f"threshold section {section!r} must be a mapping")
    allowed = {field.name for field in fields(cls)}  # type: ignore[arg-type]
    unknown = set(values) - allowed
    if unknown:
        raise SchemaError(f"unknown {section} thresholds: {sorted(unknown)}")
    try:
        if cls is GateThresholds and "required_sensitivity_analyses" in values:
            values = dict(values)
            values["required_sensitivity_analyses"] = tuple(
                values["required_sensitivity_analyses"]
            )
        return cls(**values)
    except (TypeError, ValueError) as exc:
        raise SchemaError(f"invalid {section} thresholds: {exc}") from exc


def load_thresholds(
    path: str | Path | None,
) -> tuple[LabelThresholds, EndpointThresholds, GateThresholds, dict[str, Any]]:
    """Load only declared threshold keys; defaults are explicit in returned data."""

    raw: dict[str, Any] = {}
    if path is not None:
        loaded = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        if loaded is not None and not isinstance(loaded, dict):
            raise SchemaError("threshold configuration must be a mapping")
        raw = loaded or {}
    allowed_sections = {"labeling", "endpoints", "gate", "controls", "statistics"}
    unknown_sections = set(raw) - allowed_sections
    if unknown_sections:
        raise SchemaError(f"unknown threshold sections: {sorted(unknown_sections)}")
    labels = _construct("labeling", LabelThresholds, raw.get("labeling"))
    endpoints = _construct("endpoints", EndpointThresholds, raw.get("endpoints"))
    gate = _construct("gate", GateThresholds, raw.get("gate"))

    controls = raw.get("controls") or {}
    statistics = raw.get("statistics") or {}
    if not isinstance(controls, dict) or not isinstance(statistics, dict):
        raise SchemaError("controls and statistics sections must be mappings")
    control_defaults = {
        "min_sensitivity": 0.95,
        "max_false_positive_bp_rate": 0.001,
    }
    statistics_defaults = {
        "bootstrap_iterations": 50_000,
        "permutation_iterations": 100_000,
        "seed": 20_260_729,
    }
    unknown_controls = set(controls) - set(control_defaults)
    unknown_statistics = set(statistics) - set(statistics_defaults)
    if unknown_controls:
        raise SchemaError(f"unknown controls thresholds: {sorted(unknown_controls)}")
    if unknown_statistics:
        raise SchemaError(
            f"unknown statistics thresholds: {sorted(unknown_statistics)}"
        )
    control_defaults.update(controls)
    statistics_defaults.update(statistics)
    if not 0 <= control_defaults["min_sensitivity"] <= 1:
        raise SchemaError("controls.min_sensitivity must be in [0, 1]")
    if not 0 <= control_defaults["max_false_positive_bp_rate"] <= 1:
        raise SchemaError(
            "controls.max_false_positive_bp_rate must be in [0, 1]"
        )
    if (
        not isinstance(statistics_defaults["bootstrap_iterations"], int)
        or statistics_defaults["bootstrap_iterations"] < 100
    ):
        raise SchemaError("statistics.bootstrap_iterations must be an integer >= 100")
    if (
        not isinstance(statistics_defaults["permutation_iterations"], int)
        or statistics_defaults["permutation_iterations"] < 100
    ):
        raise SchemaError(
            "statistics.permutation_iterations must be an integer >= 100"
        )
    if not isinstance(statistics_defaults["seed"], int):
        raise SchemaError("statistics.seed must be an integer")
    effective = {
        "labeling": asdict(labels),
        "endpoints": asdict(endpoints),
        "gate": asdict(gate),
        "controls": control_defaults,
        "statistics": statistics_defaults,
    }
    return labels, endpoints, gate, effective


def _load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _validate(instance: Any, schema: dict[str, Any], label: str) -> None:
    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(instance), key=lambda error: list(error.path))
    if errors:
        details = "; ".join(
            f"{'.'.join(map(str, error.path)) or '<root>'}: {error.message}"
            for error in errors
        )
        raise ValidationError(f"{label} failed schema validation: {details}")


def _resolve_root(config_path: Path) -> Path:
    config_path = config_path.resolve()
    if config_path.parent.name == "config":
        return config_path.parent.parent
    return config_path.parent


def load_and_validate(config_path: str | Path) -> ValidatedInputs:
    """Validate workflow configuration and its public sample manifest."""

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
        raise ValidationError(
            f"sample manifest columns differ; missing={missing}, extra={extra}"
        )
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
            placeholder = "0" * 32
            is_example = manifest_path.name.endswith(".example.tsv")
            if row["fastq1_md5"] != placeholder or not is_example:
                raise ValidationError(
                    f"{row['sample_id']}: mate checksums must differ"
                )

    return ValidatedInputs(
        config_path=config_path,
        root=root,
        config=config,
        manifest_path=manifest_path.resolve(),
        samples=rows,
    )
