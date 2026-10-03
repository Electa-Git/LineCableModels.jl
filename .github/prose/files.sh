#!/usr/bin/env bash
# List the maintained inputs and the pages created by Documenter.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
reports=${1:?Provide the report directory}
mkdir -p "$reports"
source_files="$reports/source-files.txt"
: > "$source_files"
while IFS= read -r -d '' file; do
    case "$file" in
        LICENSE|*/Manifest.toml|Manifest.toml|.github/prose/package-lock.json) continue ;;
        .github/prose/fixtures/*|.github/prose/styles/config/*|.github/prose/words.txt) continue ;;
        test/fixtures/reference/*/capture-*/*) continue ;;
        test/fixtures/reference/*)
            [[ "$file" == *.md || "$file" == *.jl ]] || continue ;;
    esac
    case "$file" in
        *.md|*.markdown|*.txt|*.jl|*.js|*.css|*.pro|*.sh|*.yml|*.yaml|*.json|*.toml|.gitignore)
            test -f "$file"
            printf '%s\n' "$file" >> "$source_files" ;;
    esac
done < <(git ls-files -z)
test -s "$source_files"
test -s docs/build/index.html
test -d docs/src/tutorials
find docs/src -type f -name '*.md' -print > "$reports/generated-files.txt"
grep -q '^docs/src/tutorials/.*\.md$' "$reports/generated-files.txt"
cat "$source_files" "$reports/generated-files.txt" | sort -u > "$reports/cspell-files.txt"
# Rule patterns and deliberate misspellings are data, not spelling inputs.
sed -i -e '\|^\.github/prose/styles/|d' -e '\|^\.github/prose/cspell.json$|d' "$reports/cspell-files.txt"
find docs/build -type f -name '*.html' ! -path 'docs/build/bibliography/*' ! -path 'docs/build/bibliography.html' -print > "$reports/html-files.txt"
test -s "$reports/html-files.txt"
cat "$source_files" "$reports/generated-files.txt" "$reports/html-files.txt" | sort -u > "$reports/vale-files.txt"
