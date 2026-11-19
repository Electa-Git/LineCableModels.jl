#!/usr/bin/env bash

set -Eeuo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly PLAYGROUND_DIR="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
readonly SERVER_PORT="${LCM_PRESENTATION_TEST_PORT:-18094}"
readonly DEBUG_PORT="${LCM_PRESENTATION_DEBUG_PORT:-19334}"
readonly TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/lcm-presentation.XXXXXXXX")"
server_pid=""
browser_pid=""

cleanup() {
    status=$?
    if [[ -n "${browser_pid}" ]]; then
        kill "${browser_pid}" 2>/dev/null || true
        wait "${browser_pid}" 2>/dev/null || true
    fi
    if [[ -n "${server_pid}" ]]; then
        kill "${server_pid}" 2>/dev/null || true
        wait "${server_pid}" 2>/dev/null || true
    fi
    if [[ "${status}" -ne 0 && -f "${TEST_DIR}/server.log" ]]; then
        echo "Last lines of the server's log:" >&2
        tail -n 40 "${TEST_DIR}/server.log" >&2
    fi
    rm -rf -- "${TEST_DIR}"
}
trap cleanup EXIT HUP INT TERM

# The server renders the project, then serves it. Its log goes through a pipe, which passes each line on at once.
julia --startup-file=no --project="${PLAYGROUND_DIR}" -m Playground --port "${SERVER_PORT}" \
    > >(cat >"${TEST_DIR}/server.log") 2>&1 &
server_pid=$!

for _ in {1..1200}; do
    if curl -fsS "http://127.0.0.1:${SERVER_PORT}/presentations/specimen/" >/dev/null 2>&1; then
        break
    fi
    kill -0 "${server_pid}" 2>/dev/null || { echo "the server ended before it listened" >&2; exit 1; }
    sleep 0.5
done
curl -fsS "http://127.0.0.1:${SERVER_PORT}/presentations/specimen/" >/dev/null

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

# SwiftShader gives the headless browser the WebGL of the plots.
"${browser}" --headless=new --use-angle=swiftshader --enable-unsafe-swiftshader --no-first-run \
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

LCM_PRESENTATION_ARTIFACTS="${TEST_DIR}" node "${SCRIPT_DIR}/presentation_browser.mjs" \
    "http://127.0.0.1:${SERVER_PORT}" "http://127.0.0.1:${DEBUG_PORT}" "${server_pid}"

# Each PDF has one page per slide and, in place of each cell's frame, the first line of the cell's source.
mapfile -t first_lines < <(node -e '
    for (const cell of JSON.parse(require("node:fs").readFileSync(process.argv[1], "utf8")))
        console.log(cell.source.trim().split("\n")[0]);
' "${PLAYGROUND_DIR}/_site/presentations/specimen/cells.json")
[[ "${#first_lines[@]}" == "5" ]]
for print_path in "${TEST_DIR}/direct-print.pdf" "${TEST_DIR}/preview-print.pdf"; do
    [[ "$(pdfinfo "${print_path}" | awk '/^Pages:/ {print $2}')" == "9" ]]
    text="$(pdftotext -layout "${print_path}" -)"
    for line in "${first_lines[@]}"; do
        grep -qF -- "${line}" <<<"${text}" || {
            echo "${print_path##*/} lacks the placeholder of the cell '${line}'" >&2
            exit 1
        }
    done
done

echo "Presentation integration contract passed"
