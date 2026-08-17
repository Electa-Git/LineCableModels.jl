from __future__ import annotations

import shutil
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

from .config import DEFAULT_SOURCES, read_yaml
from .util import sha256_file, utc_now, write_json


def source_catalog(path: str | Path | None = None) -> list[dict[str, Any]]:
    catalog_path = Path(path).resolve() if path else DEFAULT_SOURCES
    data = read_yaml(catalog_path)
    sources = data.get("sources", [])
    if not isinstance(sources, list):
        raise ValueError(f"Expected sources list in {catalog_path}")
    return sources


def harvest(
    root: Path,
    catalog_path: str | Path | None = None,
    *,
    download: bool = True,
    extract: bool = True,
) -> dict[str, Any]:
    root = root.resolve()
    download_root = root / "donors" / "_downloads"
    download_root.mkdir(parents=True, exist_ok=True)
    report: dict[str, Any] = {"created_at": utc_now(), "root": str(root), "sources": []}

    for entry in source_catalog(catalog_path):
        row = dict(entry)
        destination = download_root / str(entry["filename"])
        if not destination.exists():
            if not download:
                row["status"] = "missing"
                report["sources"].append(row)
                continue
            temporary = destination.with_suffix(destination.suffix + ".part")
            with urllib.request.urlopen(str(entry["url"])) as response, temporary.open("wb") as stream:
                shutil.copyfileobj(response, stream)
            temporary.replace(destination)

        actual = sha256_file(destination)
        row["path"] = str(destination)
        row["actual_sha256"] = actual
        if actual.casefold() != str(entry["sha256"]).casefold():
            row["status"] = "hash_mismatch"
            report["sources"].append(row)
            continue

        row["status"] = "verified"
        target_text = entry.get("extract_to")
        if extract and entry.get("kind") == "archive" and target_text:
            target = root / str(target_text)
            target.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(destination) as archive:
                archive.extractall(target)
            row["extracted_to"] = str(target)
        report["sources"].append(row)

    manifest = root / "donors" / "source_manifest.json"
    write_json(manifest, report)
    return report
