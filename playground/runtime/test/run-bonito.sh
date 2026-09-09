#!/usr/bin/env bash
set -Eeuo pipefail
readonly TEST_SOURCE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly RUNTIME_DIR="$(cd -- "$TEST_SOURCE/.." && pwd)"
readonly SERVER_PORT="${LCM_RUNTIME_TEST_PORT:-18116}"
readonly DEBUG_PORT="${LCM_RUNTIME_DEBUG_PORT:-19356}"
case "${LCM_RUNTIME_BROWSER_SUITE:-bonito}" in
    bonito) fixture=bonito_fixture.jl; browser_test=bonito_browser.mjs ;;
    scientific) fixture=scientific_ui_fixture.jl; browser_test=scientific_ui_browser.mjs ;;
    *) echo "Expected LCM_RUNTIME_BROWSER_SUITE=bonito or scientific" >&2; exit 2 ;;
esac
readonly TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/lcm-runtime-bonito.XXXXXXXX")"
server_pid=""
browser_pid=""
mkfifo "$TEST_DIR/shutdown"
exec 3<>"$TEST_DIR/shutdown"
cleanup() {
    result=$?
    printf 'stop\n' >&3
    for pid in "$browser_pid" "$server_pid"; do
        if [[ -n "$pid" ]]; then
            if [[ "$pid" == "$browser_pid" ]]; then kill "$pid" 2>/dev/null || true; fi
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
    exec 3>&-
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
node --check "$TEST_SOURCE/$browser_test"
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
julia --startup-file=no --threads=2 --project="$RUNTIME_DIR" \
    "$TEST_SOURCE/$fixture" "$SERVER_PORT" "$TEST_DIR" <&3 >"$TEST_DIR/server.log" 2>&1 &
server_pid=$!
wait_for "http://127.0.0.1:$SERVER_PORT/health"
browser="${LCM_BROWSER:-google-chrome}"
"$browser" --headless=new --disable-gpu --no-first-run --no-default-browser-check \
    --remote-debugging-address=127.0.0.1 --remote-debugging-port="$DEBUG_PORT" \
    --user-data-dir="$TEST_DIR/chrome" about:blank >"$TEST_DIR/chrome.log" 2>&1 &
browser_pid=$!
wait_for "http://127.0.0.1:$DEBUG_PORT/json"
node "$TEST_SOURCE/$browser_test" "http://127.0.0.1:$SERVER_PORT" "http://127.0.0.1:$DEBUG_PORT" "$TEST_DIR"
