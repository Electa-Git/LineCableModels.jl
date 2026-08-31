#!/usr/bin/env bash

set -Eeuo pipefail
IFS=$'\n\t'

readonly QUARTO_VERSION="1.9.38"
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly DATA_HOME="${XDG_DATA_HOME:-${HOME}/.local/share}"
readonly BIN_HOME="${XDG_BIN_HOME:-${HOME}/.local/bin}"
readonly QUARTO_ROOT="${DATA_HOME}/linecablemodels-playground/quarto"
readonly QUARTO_DIRECTORY="${QUARTO_ROOT}/${QUARTO_VERSION}"
readonly QUARTO_CURRENT="${QUARTO_ROOT}/current"
readonly QUARTO_EXECUTABLE="${QUARTO_DIRECTORY}/bin/quarto"

die() {
    printf 'bootstrap: %s\n' "$*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

[[ "$(uname -s)" == "Linux" ]] || die "this bootstrap currently supports Linux only"

case "$(uname -m)" in
    x86_64 | amd64)
        readonly QUARTO_ARCH="amd64"
        readonly QUARTO_SHA256="ea8c897368791ad9f200010c087ea3111b2e556b12a960487dd4e216902aa102"
        ;;
    aarch64 | arm64)
        readonly QUARTO_ARCH="arm64"
        readonly QUARTO_SHA256="75fbc5c1121ffe65e564e9d24711db2ad8f617f9552f5dc7d8a06307d72dde38"
        ;;
    *)
        die "unsupported Linux architecture: $(uname -m)"
        ;;
esac

require_command curl
require_command grep
require_command julia
require_command sha256sum
require_command tar

julia --startup-file=no -e \
    'VERSION >= v"1.12" || error("Julia 1.12 or newer is required; found $(VERSION)")'

if [[ ! -x "${QUARTO_EXECUTABLE}" ]]; then
    readonly ARCHIVE_NAME="quarto-${QUARTO_VERSION}-linux-${QUARTO_ARCH}.tar.gz"
    readonly DOWNLOAD_URL="https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/${ARCHIVE_NAME}"
    temporary_directory="$(mktemp -d "${TMPDIR:-/tmp}/linecablemodels-playground.XXXXXXXX")"
    trap 'rm -rf -- "${temporary_directory}"' EXIT
    archive_path="${temporary_directory}/${ARCHIVE_NAME}"
    extract_directory="${temporary_directory}/extract"

    printf 'Downloading Quarto %s for Linux %s...\n' "${QUARTO_VERSION}" "${QUARTO_ARCH}"
    curl --fail --location --proto '=https' --tlsv1.2 \
        --output "${archive_path}" "${DOWNLOAD_URL}"
    printf '%s  %s\n' "${QUARTO_SHA256}" "${archive_path}" | sha256sum --check --status \
        || die "Quarto archive checksum verification failed"

    mkdir -p -- "${extract_directory}" "${QUARTO_ROOT}"
    tar -xzf "${archive_path}" -C "${extract_directory}"
    extracted_directory="${extract_directory}/quarto-${QUARTO_VERSION}"
    [[ -x "${extracted_directory}/bin/quarto" ]] \
        || die "downloaded archive did not contain the expected Quarto executable"
    [[ ! -e "${QUARTO_DIRECTORY}" ]] \
        || die "incomplete Quarto installation already exists at ${QUARTO_DIRECTORY}"
    mv -- "${extracted_directory}" "${QUARTO_DIRECTORY}"
else
    printf 'Quarto %s is already installed.\n' "${QUARTO_VERSION}"
fi

if [[ -e "${QUARTO_CURRENT}" && ! -L "${QUARTO_CURRENT}" ]]; then
    die "refusing to replace non-symlink path: ${QUARTO_CURRENT}"
fi
ln -sfn -- "${QUARTO_DIRECTORY}" "${QUARTO_CURRENT}"

printf 'Instantiating the playground Julia environment...\n'
julia --startup-file=no --project="${SCRIPT_DIR}" -e \
    'using Pkg; Pkg.instantiate()'

julia_bin_directory="$(julia --startup-file=no -e 'print(joinpath(first(DEPOT_PATH), "bin"))')"
julia_command="${julia_bin_directory}/linecablemodels"
if [[ -x "${julia_command}" ]] \
    && grep -Fq "export JULIA_LOAD_PATH=${SCRIPT_DIR}" "${julia_command}"; then
    printf 'The linecablemodels application is already registered.\n'
else
    printf 'Registering the linecablemodels application...\n'
    env QUARTO_PATH="${QUARTO_EXECUTABLE}" \
        julia --startup-file=no "${SCRIPT_DIR}/install.jl"
fi
[[ -x "${julia_command}" ]] || die "Pkg did not create ${julia_command}"

mkdir -p -- "${BIN_HOME}"
public_command="${BIN_HOME}/linecablemodels"
if [[ -e "${public_command}" || -L "${public_command}" ]]; then
    if [[ ! -L "${public_command}" ]] \
        || [[ "$(readlink -- "${public_command}")" != "${julia_command}" ]]; then
        die "refusing to replace existing command: ${public_command}"
    fi
else
    ln -s -- "${julia_command}" "${public_command}"
fi

printf 'Rendering the initial Quarto site...\n'
env QUARTO_PATH="${QUARTO_EXECUTABLE}" \
    "${julia_command}" playground build --quiet

printf '\nBootstrap complete.\n'
printf '  Quarto: %s\n' "${QUARTO_EXECUTABLE}"
printf '  Command: %s\n' "${public_command}"
printf '  Start: linecablemodels playground start\n'

case ":${PATH}:" in
    *":${BIN_HOME}:"*) ;;
    *)
        printf '\n%s is not on PATH. Add this line to your shell profile:\n' "${BIN_HOME}"
        printf '  export PATH="%s:$PATH"\n' "${BIN_HOME}"
        ;;
esac
