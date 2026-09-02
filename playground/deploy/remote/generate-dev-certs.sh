#!/usr/bin/env bash

set -Eeuo pipefail
IFS=$'\n\t'
umask 077

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly OUTPUT_DIRECTORY="${1:-${SCRIPT_DIR}/certs}"
readonly SERVER_NAME="${NATS_TLS_HOST:-nats}"
readonly ARTIFACT_SERVER_NAME="${LCM_S3_TLS_HOST:-minio}"

command -v openssl >/dev/null 2>&1 || {
    printf 'generate-dev-certs: openssl was not found\n' >&2
    exit 1
}
[[ ! -e "${OUTPUT_DIRECTORY}" ]] || {
    printf 'generate-dev-certs: refusing to replace %s\n' "${OUTPUT_DIRECTORY}" >&2
    exit 1
}

temporary_directory="$(mktemp -d "${TMPDIR:-/tmp}/lcm-nats-certs.XXXXXXXX")"
trap 'rm -rf -- "${temporary_directory}"' EXIT

openssl req -x509 -newkey ec -pkeyopt ec_paramgen_curve:P-256 -nodes \
    -keyout "${temporary_directory}/ca-key.pem" \
    -out "${temporary_directory}/ca.pem" \
    -subj "/CN=LineCableModels playground development CA" \
    -days 30 >/dev/null 2>&1

issue_certificate() {
    local name="$1"
    local usage="$2"
    local extension_file="${temporary_directory}/${name}.ext"
    printf 'extendedKeyUsage=%s\n' "${usage}" >"${extension_file}"
    if [[ "${name}" == server ]]; then
        printf 'subjectAltName=DNS:%s,DNS:localhost,IP:127.0.0.1\n' \
            "${SERVER_NAME}" >>"${extension_file}"
    elif [[ "${name}" == artifact-server ]]; then
        printf 'subjectAltName=DNS:%s,DNS:localhost,IP:127.0.0.1\n' \
            "${ARTIFACT_SERVER_NAME}" >>"${extension_file}"
    fi
    openssl req -newkey ec -pkeyopt ec_paramgen_curve:P-256 -nodes \
        -keyout "${temporary_directory}/${name}-key.pem" \
        -out "${temporary_directory}/${name}.csr" \
        -subj "/O=LineCableModels/CN=${name}" >/dev/null 2>&1
    openssl x509 -req \
        -in "${temporary_directory}/${name}.csr" \
        -CA "${temporary_directory}/ca.pem" \
        -CAkey "${temporary_directory}/ca-key.pem" \
        -CAcreateserial \
        -out "${temporary_directory}/${name}-cert.pem" \
        -days 30 \
        -extfile "${extension_file}" >/dev/null 2>&1
}

issue_certificate server serverAuth
issue_certificate artifact-server serverAuth
issue_certificate publisher clientAuth
issue_certificate worker clientAuth
issue_certificate administrator clientAuth

mkdir -p -- "${OUTPUT_DIRECTORY}"
for path in \
    ca.pem \
    server-cert.pem server-key.pem \
    artifact-server-cert.pem artifact-server-key.pem \
    publisher-cert.pem publisher-key.pem \
    worker-cert.pem worker-key.pem \
    administrator-cert.pem administrator-key.pem; do
    mv -- "${temporary_directory}/${path}" "${OUTPUT_DIRECTORY}/${path}"
done

printf 'Development certificates written to %s\n' "${OUTPUT_DIRECTORY}"
printf 'They expire in 30 days and must not be committed.\n'
