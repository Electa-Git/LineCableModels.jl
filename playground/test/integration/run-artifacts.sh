#!/usr/bin/env bash

set -Eeuo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly PLAYGROUND_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
readonly CONTAINER_RUNTIME="${CONTAINER_RUNTIME:-docker}"
readonly MINIO_IMAGE="quay.io/minio/minio:RELEASE.2025-09-07T16-13-09Z@sha256:14cea493d9a34af32f524e538b8346cf79f3321eff8e708c1e2960462bd8936e"
readonly MC_IMAGE="quay.io/minio/mc:RELEASE.2025-08-13T08-35-41Z@sha256:a7fe349ef4bd8521fb8497f55c6042871b2ae640607cf99d9bede5e9bdf11727"
readonly CONTAINER_NAME="lcm-minio-integration-${$}"
readonly ROOT_ACCESS_KEY="lcm-integration-root"
readonly ROOT_SECRET_KEY="lcm-integration-root-secret"
readonly WORKER_ACCESS_KEY="lcm-integration-worker"
readonly WORKER_SECRET_KEY="lcm-integration-worker-secret"
readonly PUBLISHER_ACCESS_KEY="lcm-integration-publisher"
readonly PUBLISHER_SECRET_KEY="lcm-integration-publisher-secret"
readonly BUCKET="lcm-artifacts"

tls_mode=false
case "${1:-}" in
    "") ;;
    --tls) tls_mode=true ;;
    *)
        printf 'usage: %s [--tls]\n' "${0##*/}" >&2
        exit 2
        ;;
esac

started=false
cli_created=false
certificate_parent=""
cleanup() {
    local result=$?
    trap - EXIT HUP INT TERM
    if [[ "${cli_created}" == true ]]; then
        "${CONTAINER_RUNTIME}" rm -f "${CONTAINER_NAME}-cli" >/dev/null 2>&1 || {
            printf 'artifact integration: owned CLI cleanup failed: %s\n' "${CONTAINER_NAME}-cli" >&2
            result=1
        }
    fi
    if [[ "${started}" == true ]]; then
        "${CONTAINER_RUNTIME}" rm -f "${CONTAINER_NAME}" >/dev/null 2>&1 || {
            printf 'artifact integration: owned storage cleanup failed: %s\n' "${CONTAINER_NAME}" >&2
            result=1
        }
    fi
    if [[ -n "${certificate_parent}" ]]; then
        rm -rf -- "${certificate_parent}" || result=1
    fi
    exit "${result}"
}
trap cleanup EXIT
trap 'exit 130' HUP INT TERM

command -v "${CONTAINER_RUNTIME}" >/dev/null 2>&1 || {
    printf 'artifact integration: container runtime not found: %s\n' "${CONTAINER_RUNTIME}" >&2
    exit 1
}

container_arguments=(
    run --detach --rm
    --name "${CONTAINER_NAME}" \
    --env "MINIO_ROOT_USER=${ROOT_ACCESS_KEY}" \
    --env "MINIO_ROOT_PASSWORD=${ROOT_SECRET_KEY}" \
    --publish 127.0.0.1::9000 \
    --tmpfs /data:rw,noexec,nosuid,size=256m
)
scheme=http
minio_command=(server /data --console-address :9001)
mc_tls_arguments=()
julia_tls_environment=()
if [[ "${tls_mode}" == true ]]; then
    scheme=https
    certificate_parent="$(mktemp -d "${TMPDIR:-/tmp}/lcm-minio-tls.XXXXXXXX")"
    certificate_directory="${certificate_parent}/certs"
    "${PLAYGROUND_ROOT}/deploy/remote/generate-dev-certs.sh" \
        "${certificate_directory}" >/dev/null
    container_arguments+=(
        --volume "${certificate_directory}/artifact-server-cert.pem:/etc/minio/certs/public.crt:ro,Z"
        --volume "${certificate_directory}/artifact-server-key.pem:/etc/minio/certs/private.key:ro,Z"
    )
    minio_command=(server --certs-dir /etc/minio/certs /data --console-address :9001)
    mc_tls_arguments=(
        --volume "${certificate_directory}/ca.pem:/tmp/.mc/certs/CAs/lcm-ca.pem:ro,Z"
        --env HOME=/tmp
    )
    julia_tls_environment=("JULIA_SSL_CA_ROOTS_PATH=${certificate_directory}/ca.pem")
fi

"${CONTAINER_RUNTIME}" "${container_arguments[@]}" \
    "${MINIO_IMAGE}" "${minio_command[@]}" >/dev/null
started=true

published="$("${CONTAINER_RUNTIME}" port "${CONTAINER_NAME}" 9000/tcp)"
port="${published##*:}"
[[ "${port}" =~ ^[0-9]+$ ]] || {
    printf 'artifact integration: could not determine mapped MinIO port\n' >&2
    exit 1
}
endpoint="${scheme}://127.0.0.1:${port}"

ready=false
for _ in $(seq 1 40); do
    curl_arguments=(--fail --silent)
    [[ "${tls_mode}" == false ]] || curl_arguments+=(--cacert "${certificate_directory}/ca.pem")
    if curl "${curl_arguments[@]}" "${endpoint}/minio/health/ready" >/dev/null; then
        ready=true
        break
    fi
    sleep 0.25
done
[[ "${ready}" == true ]] || {
    printf 'artifact integration: MinIO did not become ready\n' >&2
    exit 1
}

mc() {
    # Retain an exact owned identity while attach is in flight. A shell timeout
    # must not leave an untracked CLI container behind.
    "${CONTAINER_RUNTIME}" create --name "${CONTAINER_NAME}-cli" \
        --network "container:${CONTAINER_NAME}" \
        "${mc_tls_arguments[@]}" \
        --env "MC_HOST_local=${scheme}://${ROOT_ACCESS_KEY}:${ROOT_SECRET_KEY}@127.0.0.1:9000" \
        --volume "${SCRIPT_DIR}:/policies:ro,Z" \
        "${MC_IMAGE}" "$@" >/dev/null
    cli_created=true
    "${CONTAINER_RUNTIME}" start --attach "${CONTAINER_NAME}-cli"
    "${CONTAINER_RUNTIME}" rm "${CONTAINER_NAME}-cli" >/dev/null
    cli_created=false
}

mc mb local/${BUCKET} >/dev/null
mc admin user add local "${WORKER_ACCESS_KEY}" "${WORKER_SECRET_KEY}" >/dev/null
mc admin policy create local lcm-worker /policies/minio-worker-policy.json >/dev/null
mc admin policy attach local lcm-worker --user "${WORKER_ACCESS_KEY}" >/dev/null
mc admin user add local "${PUBLISHER_ACCESS_KEY}" "${PUBLISHER_SECRET_KEY}" >/dev/null
mc admin policy create local lcm-publisher /policies/minio-publisher-policy.json >/dev/null
mc admin policy attach local lcm-publisher --user "${PUBLISHER_ACCESS_KEY}" >/dev/null

common_environment=(
    "LCM_TEST_S3_ENDPOINT=${endpoint}"
    "LCM_TEST_S3_BUCKET=${BUCKET}"
    "LCM_TEST_S3_WORKER_ACCESS_KEY=${WORKER_ACCESS_KEY}"
    "LCM_TEST_S3_WORKER_SECRET_KEY=${WORKER_SECRET_KEY}"
    "LCM_TEST_S3_PUBLISHER_ACCESS_KEY=${PUBLISHER_ACCESS_KEY}"
    "LCM_TEST_S3_PUBLISHER_SECRET_KEY=${PUBLISHER_SECRET_KEY}"
)
common_environment+=("${julia_tls_environment[@]}")

digest="$(env "${common_environment[@]}" julia --startup-file=no \
    --project="${PLAYGROUND_ROOT}/worker" \
    "${SCRIPT_DIR}/s3_worker_put.jl")"
[[ "${digest}" =~ ^[0-9a-f]{64}$ ]] || {
    printf 'artifact integration: worker returned invalid digest: %s\n' "${digest}" >&2
    exit 1
}

env "${common_environment[@]}" julia --startup-file=no \
    --project="${PLAYGROUND_ROOT}" \
    "${SCRIPT_DIR}/s3_publisher_get.jl" "${digest}"
