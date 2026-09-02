#!/usr/bin/env bash

set -Eeuo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly PLAYGROUND_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
readonly CONTAINER_RUNTIME="${CONTAINER_RUNTIME:-docker}"
readonly CONTAINER_NAME="lcm-playground-integration-${$}"
readonly NATS_IMAGE="docker.io/library/nats:2.11-alpine@sha256:e4bf19f15fd3218814a4e3c9e0064e1334bd8aa20d5984b9f1a0afd084f8cc00"
readonly PUBLISHER_PASSWORD="publisher-integration-password"
readonly WORKER_PASSWORD="worker-integration-password"
readonly ADMIN_PASSWORD="administrator-integration-password"
readonly JETSTREAM_KEY="jetstream-integration-encryption-key"

tls_mode=false
case "${1:-}" in
    "") ;;
    --tls) tls_mode=true ;;
    *)
        printf 'usage: %s [--tls]\n' "${0##*/}" >&2
        exit 2
        ;;
esac

worker_data="$(mktemp -d "${TMPDIR:-/tmp}/lcm-worker-integration.XXXXXXXX")"
certificate_parent=""
broker_started=false
broker_paused=false

cleanup() {
    if [[ "${broker_started}" == true ]]; then
        if "${CONTAINER_RUNTIME}" inspect "${CONTAINER_NAME}" \
            --format '{{.State.Paused}}' 2>/dev/null | grep -qx true; then
            "${CONTAINER_RUNTIME}" unpause "${CONTAINER_NAME}" >/dev/null 2>&1 || true
        fi
        "${CONTAINER_RUNTIME}" rm -f "${CONTAINER_NAME}" >/dev/null 2>&1 || true
    fi
    rm -rf -- "${worker_data}"
    [[ -z "${certificate_parent}" ]] || rm -rf -- "${certificate_parent}"
}

trap cleanup EXIT HUP INT TERM

command -v "${CONTAINER_RUNTIME}" >/dev/null 2>&1 || {
    printf 'integration: container runtime not found: %s\n' "${CONTAINER_RUNTIME}" >&2
    exit 1
}
command -v julia >/dev/null 2>&1 || {
    printf 'integration: Julia was not found on PATH\n' >&2
    exit 1
}

container_arguments=(
    run --detach --rm
    --name "${CONTAINER_NAME}"
    --env "NATS_PUBLISHER_PASSWORD=${PUBLISHER_PASSWORD}"
    --env "NATS_WORKER_PASSWORD=${WORKER_PASSWORD}"
    --env "NATS_ADMIN_PASSWORD=${ADMIN_PASSWORD}"
    --publish 127.0.0.1::4222
)

scheme=nats
if [[ "${tls_mode}" == true ]]; then
    scheme=tls
    certificate_parent="$(mktemp -d "${TMPDIR:-/tmp}/lcm-tls-integration.XXXXXXXX")"
    certificate_directory="${certificate_parent}/certs"
    "${PLAYGROUND_ROOT}/deploy/remote/generate-dev-certs.sh" \
        "${certificate_directory}" >/dev/null
    container_arguments+=(
        --env "NATS_JETSTREAM_KEY=${JETSTREAM_KEY}"
        --volume "${PLAYGROUND_ROOT}/deploy/remote/nats-tls.conf:/etc/nats/nats.conf:ro,Z"
        --volume "${certificate_directory}:/etc/nats/certs:ro,Z"
        --tmpfs /data:rw,noexec,nosuid,size=256m
    )
else
    container_arguments+=(
        --volume "${PLAYGROUND_ROOT}/deploy/nats.conf:/etc/nats/nats.conf:ro,Z"
    )
fi

"${CONTAINER_RUNTIME}" "${container_arguments[@]}" \
    "${NATS_IMAGE}" -c /etc/nats/nats.conf >/dev/null
broker_started=true

published="$("${CONTAINER_RUNTIME}" port "${CONTAINER_NAME}" 4222/tcp)"
port="${published##*:}"
[[ "${port}" =~ ^[0-9]+$ ]] || {
    printf 'integration: could not determine mapped NATS port from %s\n' "${published}" >&2
    exit 1
}

publisher_url="${scheme}://publisher:${PUBLISHER_PASSWORD}@127.0.0.1:${port}"
worker_url="${scheme}://worker:${WORKER_PASSWORD}@127.0.0.1:${port}"
admin_url="${scheme}://administrator:${ADMIN_PASSWORD}@127.0.0.1:${port}"
invalid_url="${scheme}://publisher:not-the-password@127.0.0.1:${port}"

publisher_tls_environment=()
admin_tls_environment=()
worker_tls_environment=()
if [[ "${tls_mode}" == true ]]; then
    publisher_tls_environment=(
        "NATS_TEST_TLS=1"
        "NATS_TLS_CA_PATH=${certificate_directory}/ca.pem"
        "NATS_TLS_CERT_PATH=${certificate_directory}/publisher-cert.pem"
        "NATS_TLS_KEY_PATH=${certificate_directory}/publisher-key.pem"
        "NATS_TLS_SERVER_NAME=localhost"
    )
    admin_tls_environment=(
        "NATS_TLS_CA_PATH=${certificate_directory}/ca.pem"
        "NATS_TLS_CERT_PATH=${certificate_directory}/administrator-cert.pem"
        "NATS_TLS_KEY_PATH=${certificate_directory}/administrator-key.pem"
        "NATS_TLS_SERVER_NAME=localhost"
    )
    worker_tls_environment=(
        "NATS_TEST_WORKER_TLS_CA_PATH=${certificate_directory}/ca.pem"
        "NATS_TEST_WORKER_TLS_CERT_PATH=${certificate_directory}/worker-cert.pem"
        "NATS_TEST_WORKER_TLS_KEY_PATH=${certificate_directory}/worker-key.pem"
        "NATS_TEST_WORKER_TLS_SERVER_NAME=localhost"
    )
fi

initialized=false
for _ in $(seq 1 20); do
    if env "${admin_tls_environment[@]}" NATS_ADMIN_URL="${admin_url}" \
        julia --startup-file=no \
        --project="${PLAYGROUND_ROOT}" \
        -e 'using LineCableModelsPlayground; LineCableModelsPlayground.initialize_runtime()'; then
        initialized=true
        break
    fi
    sleep 0.25
done
[[ "${initialized}" == true ]] || {
    printf 'integration: NATS did not become ready\n' >&2
    exit 1
}

env "${publisher_tls_environment[@]}" "${worker_tls_environment[@]}" \
    NATS_TEST_PUBLISHER_URL="${publisher_url}" \
    NATS_TEST_WORKER_URL="${worker_url}" \
    NATS_TEST_INVALID_URL="${invalid_url}" \
    NATS_TEST_WORKER_DATA="${worker_data}" \
    NATS_TEST_CONTAINER_RUNTIME="${CONTAINER_RUNTIME}" \
    NATS_TEST_CONTAINER="${CONTAINER_NAME}" \
    julia --startup-file=no --project="${PLAYGROUND_ROOT}" \
    "${PLAYGROUND_ROOT}/test/runtests.jl"
