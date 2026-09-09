#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'
readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly PLAYGROUND_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
readonly LCM_TEST_RUNTIME="${CONTAINER_RUNTIME:-podman}"
readonly LCM_TEST_SUITE="${1:-full}"
[[ "${LCM_TEST_SUITE}" == full || "${LCM_TEST_SUITE}" == terminal || "${LCM_TEST_SUITE}" == scientific || "${LCM_TEST_SUITE}" == physical ]] || exit 2
readonly LCM_TEST_CONTAINER="lcm-runtime-control-${$}"
readonly LCM_TEST_IMAGE="docker.io/library/nats:2.11-alpine@sha256:e4bf19f15fd3218814a4e3c9e0064e1334bd8aa20d5984b9f1a0afd084f8cc00"
readonly LCM_ARTIFACT_IMAGE="quay.io/minio/minio:RELEASE.2025-09-07T16-13-09Z@sha256:14cea493d9a34af32f524e538b8346cf79f3321eff8e708c1e2960462bd8936e"
readonly LCM_ARTIFACT_CLI_IMAGE="quay.io/minio/mc:RELEASE.2025-08-13T08-35-41Z@sha256:a7fe349ef4bd8521fb8497f55c6042871b2ae640607cf99d9bede5e9bdf11727"
readonly LCM_TEST_DIRECTORY="$(mktemp -d "${TMPDIR:-/tmp}/lcm-runtime-control.XXXXXXXX")"
julia_options=()
if [[ "${LCM_RUNTIME_TRACE_COMPILE:-0}" == 1 ]]; then
    julia_options+=("--trace-compile=${LCM_TEST_DIRECTORY}/compilation.jl" --trace-compile-timing)
fi
started=false
artifacts_started=false
artifacts_init_started=false
cleanup() {
    local result=$?
    trap - EXIT HUP INT TERM
    if [[ "${artifacts_init_started}" == true ]]; then
        "${LCM_TEST_RUNTIME}" rm -f "${LCM_TEST_CONTAINER}-artifacts-init" >/dev/null 2>&1 || result=1
    fi
    if [[ "${artifacts_started}" == true ]]; then
        "${LCM_TEST_RUNTIME}" logs "${LCM_TEST_CONTAINER}-artifacts" > "${LCM_TEST_DIRECTORY}/artifacts.log" 2>&1 || true
        "${LCM_TEST_RUNTIME}" rm -f "${LCM_TEST_CONTAINER}-artifacts" >/dev/null 2>&1 || result=1
    fi
    if [[ "${started}" == true ]]; then
        "${LCM_TEST_RUNTIME}" logs "${LCM_TEST_CONTAINER}" > "${LCM_TEST_DIRECTORY}/broker.log" 2>&1 || true
        "${LCM_TEST_RUNTIME}" rm -f "${LCM_TEST_CONTAINER}" >/dev/null 2>&1 || result=1
    fi
    # Remove only test credentials; retain non-secret diagnostics on failure.
    rm -f -- "${LCM_TEST_DIRECTORY}/coordinator.password" "${LCM_TEST_DIRECTORY}/worker-a.password" "${LCM_TEST_DIRECTORY}/worker-b.password"
    rm -f -- "${LCM_TEST_DIRECTORY}/artifact-coordinator.toml" "${LCM_TEST_DIRECTORY}/artifact-worker-a.toml" "${LCM_TEST_DIRECTORY}/artifact-worker-b.toml"
    if [[ -d "${LCM_TEST_DIRECTORY}/certs" && ! -L "${LCM_TEST_DIRECTORY}/certs" ]]; then
        rm -rf -- "${LCM_TEST_DIRECTORY}/certs"
    fi
    printf 'Control transport diagnostics: %s\n' "${LCM_TEST_DIRECTORY}"
    exit "${result}"
}
trap cleanup EXIT
trap 'exit 130' HUP INT TERM
"${PLAYGROUND_ROOT}/deploy/remote/generate-dev-certs.sh" "${LCM_TEST_DIRECTORY}/certs" >/dev/null
julia --startup-file=no --project="${PLAYGROUND_ROOT}/runtime" "${SCRIPT_DIR}/broker_fixture.jl" "${LCM_TEST_DIRECTORY}"
"${LCM_TEST_RUNTIME}" run --detach --rm --pull=never --name "${LCM_TEST_CONTAINER}" \
    --label lcm.runtime.test=control-v2 --publish 127.0.0.1::4222 \
    --env LCM_TEST_COORDINATOR_PASSWORD=coordinator-fixture-password \
    --env LCM_TEST_A_PASSWORD=worker-a-fixture-password \
    --env LCM_TEST_B_PASSWORD=worker-b-fixture-password \
    --env NATS_PUBLISHER_PASSWORD=publisher-fixture-password \
    --env NATS_WORKER_PASSWORD=legacy-worker-fixture-password \
    --env NATS_ADMIN_PASSWORD=administrator-fixture-password \
    --env NATS_JETSTREAM_KEY=fixture-only-jetstream-encryption-key \
    --volume "${LCM_TEST_DIRECTORY}/nats.conf:/etc/nats/nats.conf:ro,Z" \
    --volume "${LCM_TEST_DIRECTORY}/certs:/etc/nats/certs:ro,Z" \
    --tmpfs /data:rw,noexec,nosuid,size=128m \
    "${LCM_TEST_IMAGE}" -c /etc/nats/nats.conf >/dev/null
started=true
published="$("${LCM_TEST_RUNTIME}" port "${LCM_TEST_CONTAINER}" 4222/tcp)"
port="${published##*:}"
[[ "${port}" =~ ^[0-9]+$ ]] || exit 1
if [[ "${LCM_TEST_SUITE}" == physical ]]; then
    julia --startup-file=no --threads=2 "${julia_options[@]}" --project="${PLAYGROUND_ROOT}/runtime" \
        "${SCRIPT_DIR}/physical_agent_tls.jl" "${LCM_TEST_DIRECTORY}" "${port}"
    exit 0
fi
if [[ "${LCM_TEST_SUITE}" == terminal ]]; then
    julia --startup-file=no --threads=2 "${julia_options[@]}" --project="${PLAYGROUND_ROOT}/runtime" \
        "${SCRIPT_DIR}/broker_terminal.jl" "${LCM_TEST_DIRECTORY}" "${port}"
    exit 0
fi
"${LCM_TEST_RUNTIME}" run --detach --rm --name "${LCM_TEST_CONTAINER}-artifacts" \
    --label lcm.runtime.test=private-artifacts --publish 127.0.0.1::9000 \
    --env MINIO_ROOT_USER=fixture-admin --env MINIO_ROOT_PASSWORD=fixture-only-artifact-root-secret \
    --volume "${LCM_TEST_DIRECTORY}/certs/artifact-server-cert.pem:/certs/public.crt:ro,z" \
    --volume "${LCM_TEST_DIRECTORY}/certs/artifact-server-key.pem:/certs/private.key:ro,z" \
    --tmpfs /data:rw,noexec,nosuid,size=128m \
    "${LCM_ARTIFACT_IMAGE}" server /data --address :9000 --console-address :9001 --certs-dir /certs >/dev/null
artifacts_started=true
"${LCM_TEST_RUNTIME}" create --name "${LCM_TEST_CONTAINER}-artifacts-init" \
    --network "container:${LCM_TEST_CONTAINER}-artifacts" --entrypoint /bin/sh \
    --volume "${LCM_TEST_DIRECTORY}:/fixture:ro,z" \
    --volume "${SCRIPT_DIR}/artifact_initialize.sh:/initialize.sh:ro,z" \
    "${LCM_ARTIFACT_CLI_IMAGE}" /initialize.sh >/dev/null
artifacts_init_started=true
"${LCM_TEST_RUNTIME}" start --attach "${LCM_TEST_CONTAINER}-artifacts-init"
"${LCM_TEST_RUNTIME}" rm "${LCM_TEST_CONTAINER}-artifacts-init" >/dev/null
artifacts_init_started=false
artifact_published="$("${LCM_TEST_RUNTIME}" port "${LCM_TEST_CONTAINER}-artifacts" 9000/tcp)"
artifact_port="${artifact_published##*:}"
[[ "${artifact_port}" =~ ^[0-9]+$ ]] || exit 1
if [[ "${LCM_TEST_SUITE}" == scientific ]]; then
    julia --startup-file=no --threads=2 "${julia_options[@]}" --project="${PLAYGROUND_ROOT}/runtime" \
        "${SCRIPT_DIR}/broker_scientific.jl" "${LCM_TEST_DIRECTORY}" "${port}" "${artifact_port}"
    # Compare in the numerical environment, never loading an engine into the
    # coordinator. Only the complete successful browser artifact is accepted.
    timeout --signal=TERM --kill-after=10s 180s \
        julia --startup-file=no --compiled-modules=existing \
        --project="${PLAYGROUND_ROOT}/worker/profiles/line-parameters" \
        "${SCRIPT_DIR}/study_result_validation.jl" "${LCM_TEST_DIRECTORY}/scientific-results.json"
else
    julia --startup-file=no "${julia_options[@]}" --project="${PLAYGROUND_ROOT}/runtime" \
        "${SCRIPT_DIR}/broker_control.jl" "${LCM_TEST_DIRECTORY}" "${port}" "${artifact_port}"
fi
