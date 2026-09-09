#!/usr/bin/env bash
set -Eeuo pipefail
readonly TEST_SOURCE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly RUNTIME_DIR="$(cd -- "$TEST_SOURCE/.." && pwd)"
readonly PLAYGROUND_DIR="$(cd -- "$RUNTIME_DIR/.." && pwd)"
readonly SERVER_PORT="${LCM_RUNTIME_TEST_PORT:-18117}"
readonly DEBUG_PORT="${LCM_RUNTIME_DEBUG_PORT:-19357}"
readonly TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/lcm-runtime-catalogue.XXXXXXXX")"
server_pid=""
browser_pid=""
cleanup() {
    result=$?
    for pid in "$server_pid" "$browser_pid"; do
        if [[ -n "$pid" ]]; then
            kill "$pid" 2>/dev/null || true
            for _ in {1..200}; do
                kill -0 "$pid" 2>/dev/null || break
                sleep 0.1
            done
            if kill -0 "$pid" 2>/dev/null; then
                echo "Test-owned process $pid exceeded shutdown deadline" >&2
                kill -KILL "$pid" 2>/dev/null || true
                result=1
            fi
            wait "$pid" 2>/dev/null || true
        fi
    done
    if [[ -f "$TEST_DIR/runtime.sqlite" ]]; then
        julia --startup-file=no --project="$RUNTIME_DIR" "$TEST_SOURCE/catalogue_cleanup.jl" "$TEST_DIR" || result=1
    fi
    if [[ "$result" != 0 ]]; then
        for log in "$TEST_DIR/server.log" "$TEST_DIR/chrome.log"; do
            [[ ! -f "$log" ]] || tail -60 "$log" >&2
        done
    fi
    echo "Runtime browser diagnostics: $TEST_DIR"
    trap - EXIT
    exit "$result"
}
trap cleanup EXIT
trap 'exit 130' INT TERM
node --check "$TEST_SOURCE/catalogue_browser.mjs"
for port in "$SERVER_PORT" "$DEBUG_PORT"; do
    if curl -s --max-time 1 "http://127.0.0.1:$port/" >/dev/null; then
        echo "Port $port is in use; set LCM_RUNTIME_TEST_PORT / LCM_RUNTIME_DEBUG_PORT." >&2
        exit 1
    fi
done
wait_for() {
    for _ in {1..600}; do
        if curl -fsS --max-time 1 "$1" >/dev/null 2>&1; then return; fi
        sleep 0.1
    done
    echo "Timed out waiting for $1" >&2
    return 1
}
"$PLAYGROUND_DIR/lcm" playground build --quiet
julia --startup-file=no --project="$RUNTIME_DIR" "$TEST_SOURCE/catalogue_config.jl" "$SERVER_PORT" "$TEST_DIR"
"$PLAYGROUND_DIR/lcm" runtime start --config "$TEST_DIR/runtime.toml" --xray >"$TEST_DIR/server.log" 2>&1 &
server_pid=$!
wait_for "http://127.0.0.1:$SERVER_PORT/health"
browser="${LCM_BROWSER:-google-chrome}"
"$browser" --headless=new --disable-gpu --no-first-run --no-default-browser-check \
    --remote-debugging-address=127.0.0.1 --remote-debugging-port="$DEBUG_PORT" \
    --user-data-dir="$TEST_DIR/chrome" about:blank >"$TEST_DIR/chrome.log" 2>&1 &
browser_pid=$!
wait_for "http://127.0.0.1:$DEBUG_PORT/json"
node "$TEST_SOURCE/catalogue_browser.mjs" "http://127.0.0.1:$SERVER_PORT" "http://127.0.0.1:$DEBUG_PORT"
