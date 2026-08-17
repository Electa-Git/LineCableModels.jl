from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
import time
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def project_namespace(path: Path) -> str:
    with path.open("rb") as stream:
        for _, node in ET.iterparse(stream, events=("start",)):
            return node.get("name") or path.stem
    return path.stem


def project_file_version(path: Path) -> str | None:
    with path.open("rb") as stream:
        for _, node in ET.iterparse(stream, events=("start",)):
            return node.get("version")
    return None


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_json(value: Any) -> str:
    rendered = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(rendered.encode("utf-8")).hexdigest()


def normalized_lcp_bytes(path: Path) -> bytes:
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("cp1252")
    lines = [line.rstrip() for line in text.replace("\r\n", "\n").replace("\r", "\n").split("\n")]
    return ("\n".join(lines).rstrip() + "\n").encode("utf-8")


def case_id(pscad_version: str, inputs: Iterable[Path]) -> str:
    digest = hashlib.sha256()
    digest.update(pscad_version.encode("utf-8"))
    digest.update(b"\0")
    for path in sorted(inputs, key=lambda item: item.name.casefold()):
        digest.update(path.name.casefold().encode("utf-8"))
        digest.update(b"\0")
        digest.update(normalized_lcp_bytes(path))
        digest.update(b"\0")
    return digest.hexdigest()


def jsonable(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        if math.isinf(value):
            return "Infinity" if value > 0 else "-Infinity"
        # PSCAD Value objects subclass float. str(value) preserves entered units.
        if type(value) is not float:
            return str(value)
        return value
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset, tuple, list, range)):
        items = [jsonable(item) for item in value]
        return sorted(items, key=str) if isinstance(value, (set, frozenset)) else items
    if isinstance(value, Path):
        return str(value)
    return str(value)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(value), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def slug(value: str, limit: int = 52) -> str:
    result = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-") or "unnamed"
    return result[:limit]


def snapshot(root: Path) -> dict[str, tuple[int, int]]:
    if not root.exists():
        return {}
    result: dict[str, tuple[int, int]] = {}
    for path in root.rglob("*"):
        if path.is_file():
            stat = path.stat()
            result[str(path.relative_to(root))] = (stat.st_mtime_ns, stat.st_size)
    return result


def flat_snapshot(root: Path) -> dict[str, tuple[int, int]]:
    """Snapshot only direct files (used for LCP logs emitted to process CWD)."""
    result: dict[str, tuple[int, int]] = {}
    for path in root.iterdir():
        if path.is_file():
            stat = path.stat()
            result[path.name] = (stat.st_mtime_ns, stat.st_size)
    return result


def wait_for_artifacts(
    root: Path,
    before: dict[str, tuple[int, int]],
    timeout_seconds: float,
    quiet_seconds: float,
    *,
    process_cwd: Path | None = None,
    cwd_before: dict[str, tuple[int, int]] | None = None,
) -> dict[str, tuple[int, int]]:
    """Wait until canonical LCP output exists and the temp tree is stable."""
    deadline = time.monotonic() + timeout_seconds
    last = snapshot(root)
    last_change = time.monotonic()
    while time.monotonic() < deadline:
        time.sleep(0.25)
        current = snapshot(root)
        if current != last:
            last = current
            last_change = time.monotonic()
        changed = {name for name, stamp in current.items() if before.get(name) != stamp}
        if process_cwd is not None and cwd_before is not None:
            current_cwd = flat_snapshot(process_cwd)
            changed_logs = [
                process_cwd / name
                for name, stamp in current_cwd.items()
                if cwd_before.get(name) != stamp and Path(name).suffix.casefold() == ".log"
            ]
            for log in changed_logs:
                try:
                    log_text = log.read_text(encoding="cp1252", errors="replace")
                except OSError:
                    continue
                terminal_markers = (
                    "Tline.exe: ERROR",
                    "EMTTL Ending!",
                    "forrtl: severe",
                )
                if any(marker.casefold() in log_text.casefold() for marker in terminal_markers):
                    return current
        has_input = any(Path(name).suffix.lower() in {".tli", ".cli"} for name in changed)
        has_result = any(Path(name).suffix.lower() in {".tlo", ".clo", ".out"} for name in changed)
        if has_input and has_result and time.monotonic() - last_change >= quiet_seconds:
            return current
    return snapshot(root)


def copy_project_tree(source: Path, destination: Path) -> Path:
    """Copy a donor directory without PSCAD compiler/temp folders."""
    if destination.exists():
        shutil.rmtree(destination)

    def ignored(_directory: str, names: list[str]) -> set[str]:
        return {
            name
            for name in names
            if re.search(r"\.(?:gf|if|ef|gfortran|intel|clang|vc)\w*$", name, re.IGNORECASE)
            or name in {".git", "__pycache__"}
        }

    shutil.copytree(source.parent, destination, ignore=ignored)
    # Sibling cases in official multi-case archives are not dependencies. They
    # would otherwise be duplicated into every benchmark record.
    for sibling in destination.glob("*.pscx"):
        if sibling.name.casefold() != source.name.casefold():
            sibling.unlink()
    for workspace in destination.glob("*.pswx"):
        workspace.unlink()
    return destination / source.name
# End of module.
# (Kept explicit to avoid accidental executable text after this function.)
