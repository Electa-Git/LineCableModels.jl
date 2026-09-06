#!/usr/bin/env bash
set -Eeuo pipefail

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly PLAYGROUND_DIR="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
readonly SERVER_PORT="${LCM_XRAY_TEST_PORT:-18106}"
readonly DEBUG_PORT="${LCM_XRAY_DEBUG_PORT:-19346}"
readonly TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/lcm-xray.XXXXXXXX")"
server_pid=""
browser_pid=""

cleanup() {
    result=$?
    if [[ "$result" != 0 ]]; then
        for log in "$TEST_DIR/server.log" "$TEST_DIR/chrome.log"; do
            [[ ! -f "$log" ]] || sed -n '1,100p' "$log" >&2
        done
        echo "X-ray test logs retained in $TEST_DIR" >&2
    fi
    for pid in "$browser_pid" "$server_pid"; do
        if [[ -n "$pid" ]]; then
            kill "$pid" 2>/dev/null || true
            wait "$pid" 2>/dev/null || true
        fi
    done
    # Preserve diagnostics on failure; remove only this owned temporary directory.
    if [[ "$result" == 0 && "$TEST_DIR" == */lcm-xray.* ]]; then
        rm -rf -- "$TEST_DIR"
    fi
}
trap cleanup EXIT
trap 'exit 130' INT TERM

cd -- "$PLAYGROUND_DIR"
for port in "$SERVER_PORT" "$DEBUG_PORT"; do
    if curl -s --max-time 1 "http://127.0.0.1:$port/" >/dev/null; then
        echo "Test port $port is already in use; choose a different LCM_XRAY_*_PORT." >&2
        exit 1
    fi
done
./lcm playground build --quiet
julia --startup-file=no --project=. test/integration/xray_fixture.jl "$SERVER_PORT" \
    >"$TEST_DIR/server.log" 2>&1 &
server_pid=$!

wait_for() {
    for _ in {1..300}; do
        if curl -fsS --max-time 1 "$1" >/dev/null 2>&1; then return; fi
        sleep 0.1
    done
    echo "Timed out waiting for $1" >&2
    return 1
}
wait_for "http://127.0.0.1:$SERVER_PORT/"

browser="${LCM_BROWSER:-}"
if [[ -z "$browser" ]]; then
    for candidate in google-chrome chromium chromium-browser google-chrome-stable; do
        if command -v "$candidate" >/dev/null 2>&1; then
            browser="$(command -v "$candidate")"
            break
        fi
    done
fi
[[ -n "$browser" ]] || { echo "Chrome/Chromium not found" >&2; exit 1; }
"$browser" --headless=new --disable-gpu --no-first-run --no-default-browser-check \
    --remote-debugging-address=127.0.0.1 --remote-debugging-port="$DEBUG_PORT" \
    --user-data-dir="$TEST_DIR/chrome" about:blank >"$TEST_DIR/chrome.log" 2>&1 &
browser_pid=$!
wait_for "http://127.0.0.1:$DEBUG_PORT/json"
node test/integration/xray_preview_browser.mjs \
    "http://127.0.0.1:$SERVER_PORT" "http://127.0.0.1:$DEBUG_PORT"
