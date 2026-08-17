from __future__ import annotations

import argparse
import fnmatch
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any

from . import __version__
from .config import DEFAULT_CONFIG, DEFAULT_SOURCES, load_config
from .dataset import consolidate
from .runner import PscadSession, inspect_source, run_source
from .sources import harvest
from .util import project_file_version, sha256_file, utc_now, write_json


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(
        prog="linecablebenchmark",
        description="Harvest PSCAD LCP benchmark inputs, matrices, fit outputs, and diagnostics.",
    )
    root.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subcommands = root.add_subparsers(dest="command", required=True)

    harvest_parser = subcommands.add_parser("harvest", help="Verify/download and extract official PSCAD donors")
    harvest_parser.add_argument("--root", type=Path, default=Path.cwd())
    harvest_parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
    harvest_parser.add_argument("--verify-only", action="store_true", help="Do not download missing files")
    harvest_parser.add_argument("--no-extract", action="store_true")
    harvest_parser.add_argument("--json", type=Path)

    catalog_parser = subcommands.add_parser("catalog", help="Inventory local donor PSCX files and hashes")
    catalog_parser.add_argument("--root", type=Path, default=Path("donors/original"))
    catalog_parser.add_argument("--json", type=Path, default=Path("donors/project_catalog.json"))

    inspect_parser = subcommands.add_parser("inspect", help="Inspect live PSCAD datamodel and plan variants")
    add_case_arguments(inspect_parser, output=False)
    inspect_parser.add_argument("--json", type=Path)

    case_parser = subcommands.add_parser("case", help="Run all configured variants for one PSCX donor")
    add_case_arguments(case_parser, output=True)

    batch_parser = subcommands.add_parser("batch", help="Run configured variants for a donor tree")
    batch_parser.add_argument("--root", type=Path)
    batch_parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    batch_parser.add_argument("--output", type=Path)
    batch_parser.add_argument("--include", action="append", help="Override YAML include glob; repeatable")
    batch_parser.add_argument("--exclude", action="append", help="Override YAML exclude glob; repeatable")
    batch_parser.add_argument("--limit-per-case", type=int)
    batch_parser.add_argument("--max-projects", type=int)
    batch_parser.add_argument("--dry-run", action="store_true")
    batch_parser.add_argument("--no-resume", action="store_true")
    batch_parser.add_argument("--json", type=Path, default=Path("runs/batch_summary.json"))

    manifest_parser = subcommands.add_parser(
        "manifest", help="Consolidate case/batch summaries into a dataset index"
    )
    manifest_parser.add_argument(
        "--summary", type=Path, action="append", required=True, help="Run summary; repeatable"
    )
    manifest_parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    manifest_parser.add_argument("--json", type=Path, default=Path("dataset_manifest.json"))
    return root


def add_case_arguments(command: argparse.ArgumentParser, *, output: bool) -> None:
    command.add_argument("--file", type=Path, required=True, help="Donor .pscx file")
    command.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    command.add_argument("--line", help="Definition/name/IID substring; default is every unique ROW definition")
    if output:
        command.add_argument("--output", type=Path)
        command.add_argument("--dry-run", action="store_true")
        command.add_argument("--limit", type=int, help="Run only the first N variants (testing/resume aid)")
        command.add_argument("--no-resume", action="store_true")
        command.add_argument("--json", type=Path)


def inventory(root: Path) -> dict[str, Any]:
    root = root.resolve()
    projects = []
    for path in sorted(root.rglob("*.pscx"), key=lambda item: str(item).casefold()):
        projects.append(
            {
                "path": str(path),
                "relative_path": str(path.relative_to(root)),
                "size": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return {"created_at": utc_now(), "root": str(root), "project_count": len(projects), "projects": projects}


def relative_matches(path: Path, root: Path, patterns: list[str]) -> bool:
    relative = path.relative_to(root).as_posix()
    return any(fnmatch.fnmatch(relative.casefold(), pattern.casefold()) for pattern in patterns)


def batch_sources(root: Path, includes: list[str], excludes: list[str]) -> list[Path]:
    candidates = sorted(root.rglob("*.pscx"), key=lambda item: str(item).casefold())
    return [
        path
        for path in candidates
        if relative_matches(path, root, includes) and not relative_matches(path, root, excludes)
    ]


def source_config(config: dict[str, Any], source: Path) -> dict[str, Any]:
    result = deepcopy(config)
    version = project_file_version(source)
    compatibility = result["pscad"].get("compatibility_version")
    if version and version.startswith("4.") and compatibility:
        result["pscad"]["version"] = str(compatibility)
    return result


def print_result(value: dict[str, Any], destination: Path | None = None) -> None:
    if destination:
        write_json(destination.resolve(), value)
        print(destination.resolve())
    else:
        import json

        print(json.dumps(value, indent=2, sort_keys=True, default=str))


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "harvest":
            result = harvest(
                args.root,
                args.sources,
                download=not args.verify_only,
                extract=not args.no_extract,
            )
            print_result(result, args.json)
            return 1 if any(row["status"] not in {"verified"} for row in result["sources"]) else 0

        if args.command == "catalog":
            result = inventory(args.root)
            print_result(result, args.json)
            return 0

        if args.command == "manifest":
            result = consolidate(args.summary, args.config)
            print_result(result, args.json)
            return 1 if result["status_counts"].get("failed", 0) else 0

        config = load_config(args.config)
        if args.command == "inspect":
            result = inspect_source(args.file, source_config(config, args.file.resolve()))
            if args.line:
                needle = args.line.casefold()
                result["definitions"] = [
                    row
                    for row in result["definitions"]
                    if needle in row["key"]["definition"].casefold()
                    or needle in row["key"]["name"].casefold()
                    or needle == str(row["key"]["iid"])
                ]
                result["total_variants"] = sum(row["variant_count"] for row in result["definitions"])
            print_result(result, args.json)
            return 0

        if args.command == "case":
            config = source_config(config, args.file.resolve())
            result = run_source(
                args.file,
                config,
                output_root=args.output,
                line_selector=args.line,
                dry_run=args.dry_run,
                limit=args.limit,
                resume=not args.no_resume,
            )
            print_result(result, args.json)
            return 1 if result.get("failed", 0) else 0

        if args.command == "batch":
            batch_config = config["batch"]
            donor_root = (args.root or Path(batch_config["root"])) .resolve()
            includes = args.include or list(batch_config.get("include", ["**/*.pscx"]))
            excludes = args.exclude if args.exclude is not None else list(batch_config.get("exclude", []))
            sources = batch_sources(donor_root, includes, excludes)
            if args.max_projects is not None:
                sources = sources[: args.max_projects]
            output = (args.output or Path(config["output"]["root"])) .resolve()
            result: dict[str, Any] = {
                "created_at": utc_now(),
                "root": str(donor_root),
                "sources": [str(path) for path in sources],
                "runs": [],
            }
            if args.dry_run:
                for source in sources:
                    one_config = source_config(config, source)
                    one = run_source(
                        source,
                        one_config,
                        output_root=output,
                        dry_run=True,
                        limit=args.limit_per_case,
                        resume=not args.no_resume,
                    )
                    result["runs"].append(one)
            else:
                for source in sources:
                    one_config = source_config(config, source)
                    try:
                        # A fresh application session per donor prevents PSCAD
                        # Automation state/memory accumulation across hundreds
                        # of solves while retaining variant-level resume.
                        one = run_source(
                            source,
                            one_config,
                            output_root=output,
                            limit=args.limit_per_case,
                            resume=not args.no_resume,
                        )
                    except Exception as exc:
                        one = {"source": str(source), "failed": 1, "error": f"{type(exc).__name__}: {exc}"}
                    result["runs"].append(one)
            result["completed_at"] = utc_now()
            result["successful"] = sum(int(row.get("successful", 0)) for row in result["runs"])
            result["solver_rejected"] = sum(
                int(row.get("solver_rejected", 0)) for row in result["runs"]
            )
            result["provisioned"] = result["successful"] + result["solver_rejected"]
            result["failed"] = sum(int(row.get("failed", 0)) for row in result["runs"])
            print_result(result, args.json)
            return 1 if result["failed"] else 0
    except KeyboardInterrupt:
        print("Interrupted", file=sys.stderr)
        return 130
    except Exception as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    return 2
# End of module.
