#!/usr/bin/env bash
# Aggregate existing owned harnesses; never operate on a configured deployment.
set -Eeuo pipefail
IFS=$'\n\t'
umask 077
readonly SOURCE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly PLAYGROUND="$(cd -- "$SOURCE/../.." && pwd)"
readonly GROUP="${1:-all}"
usage() {
    echo 'usage: bash playground/test/integration/run-runtime-platform.sh [all|unit|browser|transport|scientific|host|physical|list]'
    echo 'CONTAINER_RUNTIME=podman|docker selects the test container engine (default podman).'
    echo 'Exit: 0 selected automated gates passed; 1 failures; 2 unavailable host prerequisites.'
    echo 'A local pass is not effective-isolation or physical remote-host certification.'
    echo 'physical is opt-in (not part of all): installed digest-pinned terminal/line/power images and LCM_TEST_CADDY are required.'
}
[[ $# -le 1 ]] || { usage >&2; exit 2; }
case "$GROUP" in
    all|unit|browser|transport|scientific|host|physical) ;;
    list) usage; exit 0 ;;
    *) usage >&2; exit 2 ;;
esac
export CONTAINER_RUNTIME="${CONTAINER_RUNTIME:-podman}"
case "$CONTAINER_RUNTIME" in podman|docker) ;; *) usage >&2; exit 2 ;; esac
if [[ "$GROUP" == physical ]]; then
    export LCM_TEST_CONTAINER_RUNTIME="$CONTAINER_RUNTIME"
    for name in LCM_TEST_TERMINAL_IMAGE LCM_TEST_LINE_IMAGE LCM_TEST_POWER_IMAGE; do
        [[ "${!name:-}" =~ ^[A-Za-z0-9][A-Za-z0-9._:/-]*@sha256:[a-f0-9]{64}$ ]] || {
            echo "$name must identify an already installed immutable image." >&2; exit 2;
        }
    done
    [[ "${LCM_TEST_CADDY:-}" == /* && -x "$LCM_TEST_CADDY" ]] || {
        echo 'LCM_TEST_CADDY must name an installed absolute Caddy executable.' >&2; exit 2;
    }
fi
command -v timeout >/dev/null
command -v flock >/dev/null
# The publication build and compiler caches are shared by this checkout. Keep
# aggregate runs serial; their own multi-owner/runtime scenarios provide load.
# Lock the existing script inode read-only, without creating a stale lock file.
exec {suite_lock}<"$SOURCE/run-runtime-platform.sh"
flock --nonblock "$suite_lock" || {
    echo 'Another aggregate run owns this checkout; wait for its report before starting another group.' >&2
    exit 2
}
readonly DIRECTORY="$(mktemp -d "${TMPDIR:-/tmp}/lcm-runtime-platform.XXXXXXXX")"
readonly REPORT="$DIRECTORY/results.tsv"
printf 'gate\tstatus\texit\tseconds\tlog\n' >"$REPORT"
printf 'Runtime platform artifacts: %s\n' "$DIRECTORY"
active=""
failed=0
unavailable=0
last_code=0
cleanup() {
    local result=$?
    trap - EXIT INT TERM HUP
    if [[ -n "$active" ]]; then
        # timeout owns this child process group and forwards TERM to it. Each
        # harness remains responsible for only its own resource cleanup.
        kill -TERM "$active" 2>/dev/null || true
        wait "$active" 2>/dev/null || true
    fi
    printf 'Retained acceptance report: %s\n' "$REPORT"
    exit "$result"
}
trap cleanup EXIT
trap 'exit 130' INT TERM HUP
cd -- "$PLAYGROUND"
{
    if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        git rev-parse HEAD
        if [[ -n "$(git status --porcelain)" ]]; then echo 'worktree=dirty'; else echo 'worktree=clean'; fi
    else
        echo 'source=exported-checkout; retain its operator-provided manifest separately'
    fi
    julia --version
    if [[ "$GROUP" != physical && "$GROUP" != host ]]; then node --version; fi
    printf 'group=%s\ncontainer_runtime=%s\n' "$GROUP" "$CONTAINER_RUNTIME"
} >"$DIRECTORY/context.txt"
selected() { [[ "$GROUP" == all || "$GROUP" == "$1" ]]; }
gate() {
    local name=$1 budget=$2 kind=$3 code=0 elapsed status started=$SECONDS
    shift 3
    printf 'RUN %s\n' "$name"
    timeout --signal=TERM --kill-after=30s "${budget}s" "$@" >"$DIRECTORY/$name.log" 2>&1 &
    active=$!
    wait "$active" || code=$?
    active=""
    elapsed=$((SECONDS-started))
    status=PASS
    if [[ "$code" != 0 ]]; then
        if [[ "$kind" == prerequisite && "$code" == 2 ]]; then status=UNAVAILABLE; unavailable=$((unavailable+1));
        else status=FAIL; failed=$((failed+1)); fi
        tail -n 35 "$DIRECTORY/$name.log" >&2
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$status" "$code" "$elapsed" "$name.log" >>"$REPORT"
    printf '%s %s (%ss)\n' "$status" "$name" "$elapsed"
    last_code=$code
    # Continue independent gates, retaining every result; the summary fails closed.
}
engine_ready=true
if selected transport || selected scientific || [[ "$GROUP" == physical ]]; then
    # Resolve through the production detector, including Docker-shim rejection.
    # Fixture transport does not require executor quotas, but must stay local.
    gate transport-engine 120 test julia --startup-file=no --project=runtime \
        runtime/test/fixture_container_engine.jl "$CONTAINER_RUNTIME" "$DIRECTORY/container-engine"
    if [[ "$last_code" == 0 ]]; then
        # Every shell/Julia fixture action, including cleanup, uses this exact
        # inspected prefix and filtered environment, never an inherited remote
        # context or a second PATH lookup of the original engine name.
        export CONTAINER_RUNTIME="$DIRECTORY/container-engine"
    else
        engine_ready=false
    fi
fi
if [[ "$GROUP" == physical && "$engine_ready" == true ]]; then
    gate physical-host 120 prerequisite ./lcm runtime check-host --runtime "$LCM_TEST_CONTAINER_RUNTIME"
    if [[ "$last_code" == 0 ]]; then
        gate physical-terminal-limits 900 test julia --startup-file=no --threads=2 --project=runtime runtime/test/physical_terminal.jl
        for stop in graceful crash; do
            gate "physical-managed-$stop" 2100 test env "LCM_TEST_AGENT_STOP=$stop" bash runtime/test/run-broker.sh physical
        done
    fi
fi
if selected unit; then
    gate fixture-engine 90 test julia --startup-file=no --project=runtime runtime/test/fixture_container_engine_test.jl
    gate protocol 300 test julia --startup-file=no --project=protocol protocol/test/runtests.jl
    gate execution-core 600 test julia --startup-file=no --threads=2 --project=worker/core worker/core/test/runtests.jl
    gate runtime 900 test julia --startup-file=no --threads=2 --project=runtime runtime/test/runtests.jl
    gate publisher 900 test julia --startup-file=no --project=. test/runtests.jl
    gate worker 1200 test env LCM_TEST_POWERIMPEDANCE=1 julia --startup-file=no --project=worker worker/test/runtests.jl
    for script in runtime_client runtime_jobs runtime_terminal; do
        gate "$script" 90 test node "test/integration/$script.mjs"
    done
fi
if selected browser; then
    gate build 600 test ./lcm playground build --quiet
    gate gateway-lifecycle 600 test julia --startup-file=no --threads=2 --project=runtime runtime/test/integration.jl
    gate graceful-shutdown 600 test julia --startup-file=no --project=. test/integration/graceful_shutdown.jl
    gate catalogue-browser 600 test bash runtime/test/run-catalogue.sh
    gate bonito-browser 600 test bash runtime/test/run-bonito.sh
    gate scientific-ui-browser 600 test env LCM_RUNTIME_BROWSER_SUITE=scientific bash runtime/test/run-bonito.sh
    gate presentations 600 test bash test/integration/run-presentation.sh
    gate ribbon 600 test bash test/integration/run-ribbon.sh
    gate xray 600 test bash test/integration/run-xray.sh
fi
if selected transport && [[ "$engine_ready" == true ]]; then
    gate runtime-tls 900 test bash runtime/test/run-broker.sh full
    gate terminal-tls 900 test bash runtime/test/run-broker.sh terminal
    gate legacy-broker-tls 900 test bash test/integration/run.sh --tls
    gate legacy-artifact-tls 600 test bash test/integration/run-artifacts.sh --tls
fi
if selected scientific && [[ "$engine_ready" == true ]]; then
    gate scientific-build 600 test ./lcm playground build --quiet
    for profile in line-parameters power-flow; do
        gate "$profile-inputs" 300 test julia --startup-file=no --compiled-modules=existing \
            "--project=worker/profiles/$profile" runtime/test/study_case_validation.jl "$profile"
    done
    gate scientific-live 2700 test bash runtime/test/run-broker.sh scientific
fi
if selected host; then
    for runtime in native podman docker; do
        gate "$runtime-prerequisites" 120 prerequisite ./lcm runtime check-host --runtime "$runtime"
    done
fi
printf 'Selected gates: %s failures, %s unavailable prerequisites.\n' "$failed" "$unavailable"
echo 'Effective executor limits, physical second-computer rehearsal and cross-runtime release certification remain separate gates; see RUNTIME_PLATFORM_PROGRESS.md.'
[[ "$failed" == 0 ]] || exit 1
[[ "$unavailable" == 0 ]] || exit 2
