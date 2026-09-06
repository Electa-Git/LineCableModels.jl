#!/usr/bin/env bash

set -Eeuo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly PLAYGROUND_DIR="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
readonly SERVER_PORT="${LCM_PRESENTATION_TEST_PORT:-18094}"
readonly DEBUG_PORT="${LCM_PRESENTATION_DEBUG_PORT:-19334}"
readonly TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/lcm-presentation.XXXXXXXX")"
readonly PDF_PATH="${TEST_DIR}/specimen.pdf"
server_pid=""
browser_pid=""

cleanup() {
    if [[ -n "${browser_pid}" ]]; then
        kill "${browser_pid}" 2>/dev/null || true
        wait "${browser_pid}" 2>/dev/null || true
    fi
    if [[ -n "${server_pid}" ]]; then
        kill "${server_pid}" 2>/dev/null || true
        wait "${server_pid}" 2>/dev/null || true
    fi
    rm -rf -- "${TEST_DIR}"
}
trap cleanup EXIT HUP INT TERM

cd -- "${PLAYGROUND_DIR}"
./lcm playground build --quiet
./lcm presentation check presentations/specimen.qmd --quiet

./lcm presentation start presentations/specimen.qmd \
    --no-render --no-open --port "${SERVER_PORT}" >"${TEST_DIR}/server.log" 2>&1 &
server_pid=$!

for _ in {1..120}; do
    if curl -fsS "http://127.0.0.1:${SERVER_PORT}/presentations/specimen.html" \
        >/dev/null 2>&1; then
        break
    fi
    sleep 0.1
done
curl -fsS "http://127.0.0.1:${SERVER_PORT}/presentations/specimen.html" \
    >/dev/null 2>&1

browser="${LCM_BROWSER:-}"
if [[ -z "${browser}" ]]; then
    for candidate in google-chrome chromium chromium-browser google-chrome-stable; do
        if command -v "${candidate}" >/dev/null 2>&1; then
            browser="$(command -v "${candidate}")"
            break
        fi
    done
fi
[[ -n "${browser}" ]] || { echo "Chrome/Chromium not found" >&2; exit 1; }

"${browser}" --headless=new --disable-gpu --no-first-run \
    --no-default-browser-check --remote-debugging-address=127.0.0.1 \
    --remote-debugging-port="${DEBUG_PORT}" \
    --user-data-dir="${TEST_DIR}/chrome" about:blank \
    >"${TEST_DIR}/chrome.log" 2>&1 &
browser_pid=$!

for _ in {1..120}; do
    if curl -fsS "http://127.0.0.1:${DEBUG_PORT}/json" >/dev/null 2>&1; then
        break
    fi
    sleep 0.1
done
curl -fsS "http://127.0.0.1:${DEBUG_PORT}/json" >/dev/null 2>&1

LCM_PRESENTATION_ARTIFACTS="${TEST_DIR}" node test/integration/presentation_browser.mjs \
    "http://127.0.0.1:${SERVER_PORT}" "http://127.0.0.1:${DEBUG_PORT}"

for print_path in "${TEST_DIR}/direct-print.pdf" "${TEST_DIR}/preview-print.pdf"; do
    [[ "$(pdfinfo "${print_path}" | awk '/^Pages:/ {print $2}')" == "9" ]]
    [[ "$(pdftotext "${print_path}" - | grep -c 'Interactive view omitted')" == "2" ]]
done

./lcm presentation export presentations/specimen.qmd --pdf \
    --output "${PDF_PATH}" --quiet
[[ "$(pdfinfo "${PDF_PATH}" | awk '/^Pages:/ {print $2}')" == "9" ]]
[[ "$(pdftotext "${PDF_PATH}" - | grep -c 'Interactive view omitted')" == "2" ]]

echo "Presentation integration contract passed"
