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

publisher_only=false
dependencies_only=false
for argument in "$@"; do
    case "${argument}" in
        --publisher-only) publisher_only=true ;;
        --dependencies-only) dependencies_only=true ;;
        -h | --help)
            printf 'Usage: %s [--publisher-only] [--dependencies-only]\n' "$0"
            exit 0
            ;;
        *) die "unknown option: ${argument}" ;;
    esac
done

if [[ "${dependencies_only}" == true ]]; then
    export JULIA_PKG_PRECOMPILE_AUTO=0
fi

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
require_command readlink
require_command sha256sum
require_command setsid
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

printf 'Instantiating the shared protocol environment...\n'
julia --startup-file=no --project="${SCRIPT_DIR}/protocol" -e \
    'using Pkg; Pkg.instantiate()'

if [[ "${publisher_only}" == false ]]; then
    printf 'Instantiating the isolated scientific worker environment...\n'
    julia --startup-file=no --project="${SCRIPT_DIR}/worker" -e \
        'using Pkg; Pkg.instantiate()'
fi

if [[ "${dependencies_only}" == true ]]; then
    printf '\nDependency bootstrap complete.\n'
    printf '  Quarto: %s\n' "${QUARTO_EXECUTABLE}"
    exit 0
fi

mkdir -p -- "${BIN_HOME}"
public_command="${BIN_HOME}/lcm"
if [[ -e "${public_command}" || -L "${public_command}" ]]; then
    if [[ ! -L "${public_command}" ]]; then
        die "refusing to replace existing command: ${public_command}"
    fi
    current_target="$(readlink -- "${public_command}")"
    if [[ "${current_target}" != "${SCRIPT_DIR}/lcm" ]]; then
        die "refusing to replace existing command: ${public_command}"
    fi
else
    ln -s -- "${SCRIPT_DIR}/lcm" "${public_command}"
fi

# Remove only the legacy launchers created by earlier revisions of this
# bootstrap. Unrelated commands or regular files are never replaced.
legacy_public_command="${BIN_HOME}/linecablemodels"
if [[ -L "${legacy_public_command}" ]]; then
    legacy_target="$(readlink -- "${legacy_public_command}")"
    if [[ "${legacy_target}" == "${SCRIPT_DIR}/linecablemodels" \
        || "${legacy_target}" == "${SCRIPT_DIR}/launcher.sh" ]]; then
        rm -- "${legacy_public_command}"
    fi
fi

julia_bin_directory="$(julia --startup-file=no -e 'print(joinpath(first(DEPOT_PATH), "bin"))')"
legacy_julia_command="${julia_bin_directory}/linecablemodels"
if [[ -f "${legacy_julia_command}" ]] \
    && grep -Fq "export JULIA_LOAD_PATH=${SCRIPT_DIR}" "${legacy_julia_command}"; then
    printf 'Removing the legacy linecablemodels Julia application...\n'
    julia --startup-file=no -e 'using Pkg; Pkg.Apps.rm("linecablemodels")'
fi

legacy_private_directory="${DATA_HOME}/linecablemodels-playground/bin"
legacy_private_command="${legacy_private_directory}/linecablemodels-julia"
if [[ -L "${legacy_private_command}" ]]; then
    rm -- "${legacy_private_command}"
    rmdir --ignore-fail-on-non-empty "${legacy_private_directory}"
fi

printf 'Rendering the initial Quarto site...\n'
env QUARTO_PATH="${QUARTO_EXECUTABLE}" \
    "${public_command}" playground build --quiet

printf '\nBootstrap complete.\n'
printf '  Quarto: %s\n' "${QUARTO_EXECUTABLE}"
printf '  Command: %s\n' "${public_command}"
printf '  Start: lcm playground start\n'
printf '  Worker: lcm worker start\n'
printf '  NATS: lcm nats status\n'
printf '  Containers: lcm container resolve\n'

case ":${PATH}:" in
    *":${BIN_HOME}:"*) ;;
    *)
        printf '\n%s is not on PATH. Add this line to your shell profile:\n' "${BIN_HOME}"
        printf '  export PATH="%s:$PATH"\n' "${BIN_HOME}"
        ;;
esac
