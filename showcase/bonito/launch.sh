#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repository_directory="$(cd -- "${script_directory}/../.." && pwd)"
target="${1:-/}"

cd -- "${repository_directory}"
exec julia --threads=auto --project="${script_directory}" "${script_directory}/serve.jl" --open "${target}"
