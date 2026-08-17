from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml


PACKAGE_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = PACKAGE_ROOT / "benchmark.yml"
DEFAULT_SOURCES = PACKAGE_ROOT / "official_sources.yml"


def merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    result = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = deepcopy(value)
    return result

def read_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        value = yaml.safe_load(stream) or {}
    if not isinstance(value, dict):
        raise ValueError(f"Expected a YAML mapping in {path}")
    return value


def load_config(path: str | Path | None) -> dict[str, Any]:
    defaults = read_yaml(DEFAULT_CONFIG)
    if path is None:
        result = defaults
    else:
        result = merge(defaults, read_yaml(Path(path).resolve()))
    result["_config_path"] = str(Path(path).resolve() if path else DEFAULT_CONFIG)
    return result
