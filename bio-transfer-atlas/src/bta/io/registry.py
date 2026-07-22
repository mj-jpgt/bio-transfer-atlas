"""
Dataset registry loader — reads configs/datasets.yaml and resolves local paths.
"""

from pathlib import Path
from typing import Any

import yaml


def load_registry(config_path: Path | None = None) -> dict[str, Any]:
    if config_path is None:
        config_path = Path(__file__).resolve().parents[4] / "configs" / "datasets.yaml"
    with open(config_path) as f:
        return yaml.safe_load(f)


def get_local_path(key: str, registry: dict | None = None) -> Path:
    if registry is None:
        registry = load_registry()
    entry = registry[key]
    return Path(entry["local"])


def get_url(key: str, registry: dict | None = None) -> str:
    if registry is None:
        registry = load_registry()
    return registry[key]["url"]
