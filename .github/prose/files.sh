#!/usr/bin/env bash
# List the maintained inputs.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
reports=${1:?Provide the report directory}
mkdir -p "$reports"
source_files="$reports/source-files.txt"
: > "$source_files"
while IFS= read -r -d '' file; do
    case "$file" in
        LICENSE|*/Manifest.toml|Manifest.toml|.github/prose/package-lock.json) continue ;;
        .github/prose/fixtures/*|.github/prose/styles/config/*) continue ;;
    esac
    case "$file" in
        *.md|*.markdown|*.txt|*.jl|*.js|*.css|*.sh|*.yml|*.yaml|*.json|*.toml|.gitignore)
            test -f "$file"
            printf '%s\n' "$file" >> "$source_files" ;;
    esac
done < <(git ls-files -z)
test -s "$source_files"
sort -u "$source_files" > "$reports/vale-files.txt"
# Rule patterns and deliberate misspellings are data, not spelling inputs.
sed -e '\|^\.github/prose/styles/|d' -e '\|^\.github/prose/cspell.json$|d' \
    "$reports/vale-files.txt" > "$reports/cspell-files.txt"
