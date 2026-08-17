from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from .config import load_config
from .util import sha256_file, utc_now


def case_summaries(value: dict[str, Any]) -> Iterable[dict[str, Any]]:
    runs = value.get("runs")
    if isinstance(runs, list):
        yield from (row for row in runs if isinstance(row, dict))
    else:
        yield value


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(Path.cwd().resolve()))
    except ValueError:
        return str(path)


def consolidate(summary_paths: list[Path], config_path: Path) -> dict[str, Any]:
    """Create one portable index from case and/or batch run summaries."""
    summaries: list[dict[str, Any]] = []
    records: dict[str, dict[str, Any]] = {}
    excluded: list[dict[str, Any]] = []

    for summary_path in summary_paths:
        resolved_summary = summary_path.resolve()
        value = json.loads(resolved_summary.read_text(encoding="utf-8"))
        summaries.append(
            {"path": display_path(resolved_summary), "sha256": sha256_file(resolved_summary)}
        )
        for case in case_summaries(value):
            source = case.get("source")
            for row in case.get("excluded_definitions", []):
                excluded.append({"source": source, **row})
            for row in case.get("records", []):
                record_dir = Path(row["path"]).resolve()
                key = str(record_dir).casefold()
                manifest_path = record_dir / "manifest.json"
                if not manifest_path.is_file():
                    records[key] = {
                        "path": display_path(record_dir),
                        "status": "missing_manifest",
                    }
                    continue
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                records[key] = {
                    "path": display_path(record_dir),
                    "status": manifest.get("status"),
                    "case_id": manifest.get("case_id"),
                    "specification_sha256": manifest.get("specification_sha256"),
                    "source": manifest.get("source"),
                    "source_sha256": manifest.get("source_sha256"),
                    "pscad_version": manifest.get("pscad_version"),
                    "line": manifest.get("line"),
                    "variant": manifest.get("variant"),
                    "artifact_count": len(manifest.get("artifacts", [])),
                }

    ordered = sorted(records.values(), key=lambda row: str(row["path"]).casefold())
    counts = Counter(str(row.get("status", "unknown")) for row in ordered)
    canonical_by_case: dict[str, str] = {}
    duplicate_groups: dict[str, dict[str, Any]] = {}
    for row in ordered:
        identifier = row.get("case_id")
        if not identifier:
            continue
        identifier = str(identifier)
        canonical = canonical_by_case.get(identifier)
        if canonical is None:
            canonical_by_case[identifier] = str(row["path"])
            continue
        row["effective_duplicate_of"] = canonical
        group = duplicate_groups.setdefault(
            identifier,
            {"case_id": identifier, "canonical": canonical, "aliases": []},
        )
        group["aliases"].append(str(row["path"]))
    config_resolved = config_path.resolve()
    return {
        "schema_version": 1,
        "created_at": utc_now(),
        "config": {
            "path": display_path(config_resolved),
            "sha256": sha256_file(config_resolved),
            "effective": load_config(config_resolved),
        },
        "input_summaries": summaries,
        "record_count": len(ordered),
        "unique_effective_case_count": len(canonical_by_case),
        "duplicate_effective_record_count": sum(
            len(group["aliases"]) for group in duplicate_groups.values()
        ),
        "duplicate_effective_case_groups": sorted(
            duplicate_groups.values(), key=lambda row: row["case_id"]
        ),
        "status_counts": dict(sorted(counts.items())),
        "excluded_definition_count": len(excluded),
        "excluded_definitions": excluded,
        "records": ordered,
    }
