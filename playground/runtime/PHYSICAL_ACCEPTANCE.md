# Physical Docker / Podman acceptance

This is an opt-in test on an operator-owned Linux account with working user
systemd, cgroup v2 CPU/memory/PID delegation, and an approved local engine. It
does not configure a public deployment. Do not use the local IT-managed machine
to install Docker or weaken missing limits.

Use a separate checkout/private depot and no existing `lcm-agent-worker-a.service`.
The harness refuses that unit rather than replacing a configured worker. Keep
the account free of concurrent container mutations while testing: exact before/
after inventories form part of the cleanup evidence. Unrelated resources are
never pruned. The current owned-host results live in
[RUNTIME_PLATFORM_PROGRESS.md](../RUNTIME_PLATFORM_PROGRESS.md).

## Prepare explicitly

Install Julia 1.12.7, the selected engine, OpenSSL and an approved Caddy executable.
Instantiate `playground/runtime` in the private depot. Image construction is an
operator action; admission and registration never build or pull images.

From the repository root, for example with Podman:

```bash
for profile in julia-terminal line-parameters power-flow; do
    podman build --memory=4g --memory-swap=4g \
        --cpu-period=100000 --cpu-quota=200000 \
        --target "$profile" -f playground/worker/containers/Containerfile \
        -t "localhost/lcm-$profile:acceptance" .
done
podman image inspect localhost/lcm-julia-terminal:acceptance --format '{{json .RepoDigests}}'
podman image inspect localhost/lcm-line-parameters:acceptance --format '{{json .RepoDigests}}'
podman image inspect localhost/lcm-power-flow:acceptance --format '{{json .RepoDigests}}'
```

Use `docker build`/`docker image inspect` for Docker. BuildKit builder resource
configuration is separate from executor limits; the recipe bounds Julia package
precompilation to two tasks on either engine. Select the exact manifest references
from `RepoDigests`, not mutable tags or image config IDs. Check the actual image
inventory if the builder requires an explicit local export.

Pre-pull the pinned fixture broker on the selected engine (the harness uses
`--pull=never`):

```bash
podman pull docker.io/library/nats:2.11-alpine@sha256:e4bf19f15fd3218814a4e3c9e0064e1334bd8aa20d5984b9f1a0afd084f8cc00
```

## Run the same gate on each engine

Set the three references to the corresponding engine's installed images:

```bash
export CONTAINER_RUNTIME=podman
export LCM_TEST_TERMINAL_IMAGE='localhost/lcm-julia-terminal@sha256:REPLACE_WITH_64_HEX_DIGEST'
export LCM_TEST_LINE_IMAGE='localhost/lcm-line-parameters@sha256:REPLACE_WITH_64_HEX_DIGEST'
export LCM_TEST_POWER_IMAGE='localhost/lcm-power-flow@sha256:REPLACE_WITH_64_HEX_DIGEST'
export LCM_TEST_CADDY=/absolute/private/tools/caddy
bash playground/test/integration/run-runtime-platform.sh physical
```

Repeat with `CONTAINER_RUNTIME=docker` and Docker's manifest references. The
aggregate resolves and pins its local broker-engine command; the agent separately
verifies its own local executor engine. This explicit group is not silently run
by `all`, because it starts real executors and intentionally exhausts their
bounded memory/PID budgets. Missing image references or Caddy fail before launch.

The group runs effective kernel/PTY limits, then real managed-agent TLS/Caddy
tests with graceful stop and deliberate agent death. Both owners keep isolated
REPL state while the approved line-parameter and power-flow profiles prepare and
execute through the protected job API. The coordinator loads no scientific
packages. A 202 response means acceptance, not readiness or success.

Each run reports private diagnostic directories and the aggregate `results.tsv`.
Credentials/certificates are removed by the fixture; non-secret diagnostics stay
for inspection. Exact receipt verification and the fixed systemd post-stop
recovery command own executor cleanup, including failed tests. A failed recovery
is a failure, never an excuse for a broad container prune.

This gate does not itself prove physical second-computer transport, browser
rendering, public DNS/certificates, or production artifact storage. Record those
separately. A successful image build or native transport fixture is not physical
container acceptance.

## Separate-computer consumer gate

`runtime/test/physical_consumers.jl` is a second, opt-in gate. It runs the existing
Chrome/Bonito/Reveal consumer tests on the gateway computer against already
provisioned managed agents on another computer. It neither installs those agents
nor exposes an SSH orchestration endpoint. Use a fresh private fixture deployment,
not production accounts or services.

The operator supplies:

- A private coordinator configuration for TLS NATS and private HTTPS S3 storage.
  The two machines have separate scratch roots and SQLite/process identities; do
  not mount the executor filesystem into the gateway.
- Approved, running `worker-a` (line-parameters and julia-terminal, capacity 4)
  and `worker-b` (power-flow, capacity 2), with immutable container profiles and
  4 GiB / 600 s preparation / 120 s job budgets for the scientific specimen.
  The agent engine may be Docker or Podman; its installed digests must match the
  coordinator declarations. Install the actual rendered managed service contract.
- An existing publication build and installed UI/profile environments on the
  gateway computer. The coordinator itself must not import those UI or numerical
  packages. Node and owned headless Chrome are required for this gate.
- A mode-0600 JSON array of operator-owned command arguments that stops **only**
  this fixture's power-flow service and reconciles its exact acquisition journal.
  The first argument is an absolute executable. This private test input is never
  accepted by the browser or copied into application configuration.

From the repository root on the gateway computer:

```bash
LCM_PHYSICAL_ARTIFACT_CHECK=1 julia --startup-file=no --threads=2 \
  --project=playground/runtime playground/runtime/test/physical_consumers.jl \
  /absolute/private/control.toml /absolute/private/new-diagnostics \
  /absolute/private/stop-owned-power-worker.json

julia --startup-file=no --compiled-modules=existing \
  --project=playground/worker/profiles/line-parameters \
  playground/runtime/test/study_result_validation.jl \
  /absolute/private/new-diagnostics/scientific-results.json
```

The browser measures both real calculations in both application hosts, explicit
cold/warm preparation, cancellation/retry, retained results after worker loss,
continued independent line execution, and a result large enough to require
private S3 transport. The same real terminal renderer then exercises the remote
REPL, including reconnect/restart, completion, Unicode, resize and both themes.
No fake browser responses, fabricated numerical payloads or shared result files
stand in for those transports. Only the numerical parity check consumes the
result evidence file afterwards, in a separate engine process.

Private Tailscale/SSH forwarding may provide the operator's test connectivity
where a rootless Tailscale adapter requires it. This does not change the runtime
architecture: agents initiate broker connections, executors have no network,
and there is no inbound Julia port. Verified TLS remains enabled inside that
private connection. The loopback development identity used by the browser gate
is not a public authentication deployment; the separate physical Caddy gate
tests the protected deployment boundary.

The caller must finally stop/reconcile both exact owned agents, close its owned
forwarding process, remove only its fixture services, and delete the test
credentials/certificates. Retain non-secret logs, screenshots, timings and
results, and compare remote inventory against baseline. Failure to complete
cleanup fails the rehearsal even when browser assertions passed.

The full gate uses a separate, private operator fault-driver process and waits
for its verified exit. It keeps the real presentation connected during the
companion terminal test; neither test-only process compilation nor post-run
numerical comparisons should consume a live lease's acknowledgement budget.
For a focused terminal regression, the same harness accepts
`LCM_PHYSICAL_CONSUMERS=terminal`: only `worker-a` must be online, no scientific
job or worker-stop command is issued, and no complete consumer/S3 pass is claimed.
Use another fresh diagnostic directory for every mode or engine repetition.

The terminal browser receives the shipped 120-second guarded startup bound for
Connect and explicit Restart; normal DOM interactions and reconnect keep their
short 7.5-second bound. The parent permits both startup periods plus a finite
interaction budget, and the companion presentation stays connected for that
bounded test. Failed/exited/uncertain/disconnected terminal states fail promptly.
These are fixture deadlines, not changes to production startup/lease policy.
The gate prints measured guarded startup time; do not mistake an open socket
or an accepted request for a ready Julia prompt.

The separate numerical comparison uses default Float64 relative tolerance and
a reference-slice-scaled machine-epsilon absolute floor for cancellation
residuals near zero. It includes negative controls rejecting meaningful errors,
NaN and nonzero values at zero reference scale; there is no fixed tolerance
shared between impedance and admittance units.
