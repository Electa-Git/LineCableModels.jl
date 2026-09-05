#!/usr/bin/env python3
"""Validate a collection release against Git tags and package staged benchmarks."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


COLLECTION_PATTERN = re.compile(r"^[a-z][a-z0-9_]*$")
VERSION_PATTERN = re.compile(r"^(\d+)\.(\d+)\.(\d+)$")


def run_git(repository: Path, *arguments: str) -> str:
    result = subprocess.run(
        ("git", *arguments),
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def parse_version(value: str) -> tuple[int, int, int]:
    match = VERSION_PATTERN.fullmatch(value)
    if match is None:
        raise argparse.ArgumentTypeError("version must have the form MAJOR.MINOR.PATCH")
    return tuple(int(part) for part in match.groups())


def format_version(version: tuple[int, int, int]) -> str:
    return ".".join(str(part) for part in version)


def released_versions(repository: Path, collection: str) -> list[tuple[int, int, int]]:
    prefix = f"gauntlet-{collection}-v"
    tags = run_git(repository, "tag", "--list", f"{prefix}*").splitlines()
    versions: list[tuple[int, int, int]] = []
    for tag in tags:
        suffix = tag.removeprefix(prefix)
        match = VERSION_PATTERN.fullmatch(suffix)
        if match is not None:
            versions.append(tuple(int(part) for part in match.groups()))
    return sorted(set(versions))


def allowed_successors(version: tuple[int, int, int]) -> set[tuple[int, int, int]]:
    major, minor, patch = version
    return {
        (major, minor, patch + 1),
        (major, minor + 1, 0),
        (major + 1, 0, 0),
    }


def validate_target(
    target: tuple[int, int, int],
    prior: list[tuple[int, int, int]],
) -> tuple[int, int, int] | None:
    if not prior:
        if target != (1, 0, 0):
            raise ValueError("the first release of a collection must be 1.0.0")
        return None
    latest = prior[-1]
    if target not in allowed_successors(latest):
        allowed = ", ".join(
            format_version(version) for version in sorted(allowed_successors(latest))
        )
        raise ValueError(
            f"{format_version(target)} is not a direct successor of "
            f"{format_version(latest)}; expected one of: {allowed}"
        )
    return latest


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Package one staged Gauntlet collection after validating its explicit "
            "version against collection-specific Git tags."
        )
    )
    parser.add_argument("--collection", required=True)
    parser.add_argument("--version", required=True, type=parse_version)
    parser.add_argument("--reason", required=True)
    parser.add_argument("--julia", default="julia", help="Julia executable")
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace an existing local package for this exact version",
    )
    parser.add_argument(
        "--create-tag",
        action="store_true",
        help="create the validated annotated Git tag after packaging",
    )
    return parser.parse_args()


def main() -> int:
    options = arguments()
    repository = Path(__file__).resolve().parents[2]
    collection = options.collection
    if COLLECTION_PATTERN.fullmatch(collection) is None:
        raise ValueError(
            "collection must use lowercase letters, digits, and underscores"
        )
    reason = options.reason.strip()
    if not reason:
        raise ValueError("release reason cannot be empty")
    dirty = run_git(repository, "status", "--porcelain")
    if dirty:
        raise RuntimeError(
            "the worktree is dirty; commit the definitive benchmark sources before "
            "packaging"
        )
    versions = released_versions(repository, collection)
    latest = validate_target(options.version, versions)
    version = format_version(options.version)
    commit = run_git(repository, "rev-parse", "HEAD")
    tag = f"gauntlet-{collection}-v{version}"
    command = [
        options.julia,
        "--project=test/gauntlet",
        "--startup-file=no",
        "test/gauntlet/package.jl",
        collection,
        version,
        reason,
        commit,
        str(options.force).lower(),
    ]
    print(
        f"Packaging {collection} {version}; previous release: "
        f"{format_version(latest) if latest is not None else 'none'}",
        flush=True,
    )
    subprocess.run(command, cwd=repository, check=True)
    if options.create_tag:
        run_git(
            repository,
            "tag",
            "--annotate",
            tag,
            commit,
            "--message",
            f"Gauntlet {collection} {version}: {reason}",
        )
        print(f"Created local annotated tag {tag}")
    else:
        print(f"Validated release tag: {tag}")
        print("No Git tag was created; pass --create-tag when that is intended.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, subprocess.CalledProcessError, RuntimeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(2)
