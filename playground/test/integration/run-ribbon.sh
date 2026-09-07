#!/usr/bin/env bash
set -Eeuo pipefail

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly PLAYGROUND_DIR="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
readonly SERVER_PORT="${LCM_RIBBON_TEST_PORT:-18102}"
readonly DEBUG_PORT="${LCM_RIBBON_DEBUG_PORT:-19342}"
readonly TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/lcm-ribbon.XXXXXXXX")"
server_pid=""
browser_pid=""

cleanup() {
    result=$?
    if [[ "$result" != 0 ]]; then
        for log in "$TEST_DIR/server.log" "$TEST_DIR/chrome.log"; do
            [[ ! -f "$log" ]] || sed -n '1,100p' "$log" >&2
        done
    fi
    for pid in "$browser_pid" "$server_pid"; do
        if [[ -n "$pid" ]]; then
            kill "$pid" 2>/dev/null || true
            for _ in {1..50}; do
                kill -0 "$pid" 2>/dev/null || break
                sleep 0.1
            done
            # Only processes started by this test runner are owned here.
            kill -0 "$pid" 2>/dev/null && kill -KILL "$pid" 2>/dev/null || true
            wait "$pid" 2>/dev/null || true
        fi
    done
    rm -rf -- "$TEST_DIR"
}
trap cleanup EXIT
trap 'exit 130' INT TERM

cd -- "$PLAYGROUND_DIR"
for port in "$SERVER_PORT" "$DEBUG_PORT"; do
    if curl -s --max-time 1 "http://127.0.0.1:$port/" >/dev/null; then
        echo "Test port $port is in use; choose different LCM_RIBBON_* ports." >&2
        exit 1
    fi
done
./lcm playground build --quiet
julia --startup-file=no --project=. test/integration/ribbon_fixture.jl "$SERVER_PORT" \
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

node test/integration/ribbon_theme_browser.mjs \
    "http://127.0.0.1:$SERVER_PORT" "http://127.0.0.1:$DEBUG_PORT"
echo "Shared ribbon/toolbar theme contract passed"
