# Runtime platform v1 — implementation ledger

Goal started: 2026-09-07. Contract: [RUNTIME_PLATFORM_PLAN.md](RUNTIME_PLATFORM_PLAN.md).

## Checkpoints

| Package | Status | Evidence |
|---|---|---|
| P0 baseline/configuration inventory | Complete | Unit and browser gates; pre-existing timing race recorded below |
| P1 ownership and UI-host isolation | Complete | 101 unit + 92 HTTP/WebSocket/process assertions; real Bonito browser/reconnect/restart gate passes |
| P2 public/developer catalogue | Complete | Actual CLI/catalogue/upload/teardown browser gate, 132 unit + 100 integration assertions, publisher/presentation/ribbon regressions pass |
| P3 worker control and targeting | Complete | Fresh complete TLS transport group passes: control/leases, targeted delivery, private relay, authorization, cancellation/recovery and legacy compatibility. Shared client/component conformance and registered-host browser gates pass; physical executor isolation remains P4/P5/P7 |
| P4 profiles and preparation | Complete | Real line-parameter/power-flow profiles pass on both rootless Docker and Podman, through actual Caddy/TLS with graceful stop and forced agent death. Exact cleanup, fresh readiness and kernel isolation are verified |
| P5 private Julia terminal | Complete | Shared JuliaTerminal and actual Chrome/relay/Julia path verified; both physical engines pass independent-owner REPLs, enforced memory/PID limits, flood containment and exact cleanup |
| P6 scientific consumers | Complete | Shared showcase/CableStudy cases pass native and separate-computer Docker/Podman gates. Each remote full run passes 49 scientific browser + 15 terminal browser + 29 coordinator checks, followed by 39 direct-engine comparisons and 8 roundoff-tolerance controls |
| P7 release verification | Complete | All mandatory local unit/browser/transport/offline gates, both physical engine matrices, actual Caddy authorization and both separate-computer browser/S3/terminal runs pass. Latest runtime: 2,289 assertions; final TLS terminal: 435; terminal rendering: 42. Exact local/remote teardown audits pass. No public deployment is claimed |

## Baseline

Repository: `LineCableModels-playground`, branch `release/v0.2.0-playground`,
commit `d50367bfed5174b77882cffa9116ce5b3697c6b2`. Initial worktree contained
only the approved, untracked execution plan. No branch switch or reset.

Toolchain: Julia 1.12.7, Node 24.0.0, installed Chrome and pinned Quarto.
Rootless Podman 5.8.2 uses cgroups v2 with `memory` and `pids`, **not cpu**.
The `docker` executable is a Podman compatibility shim; actual Docker
integration remains unverified. No existing containers were changed.

The first full baseline logs are retained under
`/tmp/lcm-runtime-platform.9uEfiTvj/` on this machine. This temporary location
is diagnostic evidence, not a deployment path or a substitute for repeatable
commands.

| Gate | Result |
|---|---|
| `julia --startup-file=no --project=playground/protocol playground/protocol/test/runtests.jl` | Pass: 43 assertions |
| `julia --startup-file=no --project=playground playground/test/runtests.jl` | Pass: all groups, including 433 visual-contract assertions |
| `julia --startup-file=no --project=playground/worker playground/worker/test/runtests.jl` | Pass: 70 assertions, including real engine isolation/parity |
| `bash playground/test/integration/run-presentation.sh` | Pass: published shell, fragments, math notes, geometry and PDF |
| `bash playground/test/integration/run-ribbon.sh` | Initial hover-frame race; diagnostic rerun and settled-input rerun pass |
| `bash playground/test/integration/run-xray.sh` | Pass: complete existing browser gate |

Worker baseline emitted existing precompilation/version warnings and one
duplicate-docstring warning. These did not fail its tests. Optional
`LCM_TEST_POWERIMPEDANCE=1` coverage still needs its separate release run.

### Original environment checksums (SHA-256)

| Environment file | Digest |
|---|---|
| `playground/Project.toml` | `224e9ba02c3c0b0a69c59268061245e3c5a91a3026108fe71e6e18f31ddc83aa` |
| `playground/Manifest.toml` | `f7638e0650f06a91d823b3f1b5cedf04b994cd0deb1e8f08d4393f11deaacce3` |
| `playground/worker/Project.toml` | `5487c7169e3c60d6e0502263084bd71c91416bd9711f67ecb82e9e2bd71bfeff` |
| `playground/worker/Manifest.toml` | `5e8598973ac3738b06aa8a0e9d33b18510b4df1c031a270b98289b4ee09419b5` |
| `playground/protocol/Project.toml` | `fde29719c3de889d41f4e95e70b36623127b740fdddfe96018221724fd69f090` |
| `playground/protocol/Manifest.toml` | `78156f407792c28de25c5d61e303c169d2d40339b2582b7ffba171194ead9939` |

## Frozen migration inventory

Reference old branch commit: `5bc9916a`. Port narrative/inputs/assets, not
legacy page/resource or REPL execution infrastructure.

| ICHQP2026 source section | Destination responsibility |
|---|---|
| `010_opening.jl` | Quarto opening and research context |
| `020_fundamentals.jl` | Quarto fundamentals and lightweight visual case |
| `030_cableconstants.jl` | Cable geometry/constants view |
| `040_lineparameters.jl` | Shared line-parameter view and registered operations |
| `050_applicationcase.jl` | OHL/UGC case, separate power-flow preparation |
| `060_conclusion.jl` | Quarto conclusions |

One `showcase.qmd` deck, one `CableStudy` workbench, one portable terminal.
Existing `WIDGET_ROUTES`, `WORKBENCH_ROUTES`, toolkit component render methods
and X-ray declarations are the starting conformance inventory. The inventory
must be derived from registration when new components are added.

## Requirement → verification matrix

| ID | Requirement | Acceptance evidence to add |
|---|---|---|
| RT01 | Versioned config; no implicit execution | Strict config/catalogue tests, offline startup |
| RT02 | Required dispatch hooks; cheap registration | Missing-hook/duplicate/version tests |
| RT03 | Run ownership; proxy identity; origin/CSRF | Two-principal HTTP/WebSocket denial tests |
| RT04 | Independent Bonito hosts and persistent frames | Browser identity/theme/X-ray + host-death test |
| RT05 | Bounded startup/admission/cleanup | Process inventory and repeated failure tests |
| RT06 | Durable approved registration, ephemeral health | SQLite restart + authenticated announcement tests |
| RT07 | Targeted, fenced leases/jobs | Competing agents, stale boot/grant/result tests |
| RT08 | Independent, honest preparation | Cold/warm/replacement/cache/cancel tests |
| RT09 | Bounded correlated diagnostics | Flood, redaction, disconnect-gap tests |
| RT10 | Private real Julia PTY | Completion/multiline/resize/interrupt/reconnect tests |
| RT11 | Enforced terminal limits and isolation | Actual runtime preflight, adversarial subprocess tests |
| RT12 | Shared components and both themes | Registered-host visual/interaction/X-ray matrix |
| RT13 | Restored scientific frontface and consumers | Offline catalogue + numerical parity/browser tests |
| RT14 | Compatibility and release recovery | Existing suites + aggregate fault-injection runner |

## Original external verification constraints

These describe the starting environment, not current release status. Both
engines have now passed the owned-host physical matrices recorded below; the
local IT-managed host remains unchanged and correctly fails CPU admission.

- This host lacks a delegated CPU controller. Terminal launch must fail closed
  for profiles requiring CPU quotas; do not remove that requirement to pass.
- No real Docker Engine/runner has been identified. Implement adapter tests,
  but retain actual Docker verification as pending until a suitable runner is
  available.
- No public deployment, infrastructure reconfiguration, or unrelated container
  cleanup is part of this goal.

## P1 implementation evidence

The new `runtime/` package has a separate resolved environment with no Bonito
or scientific dependency. `@required` validates application declaration/launch
hooks; immutable catalogue metadata is captured once. Strict configuration
rejects unsafe listeners, broad scratch paths and unsupported versions.

Proxy identity requires an approved actual peer, a private credential and a
principal. Mutations/WebSockets additionally check the configured Origin;
mutations require a non-simple request header. Local development also rejects
foreign Host authorities to prevent a loopback DNS-rebinding bypass. Real HTTP
and text/binary WebSocket checks pass with two distinct proxy principals.

SQLite ownership/admission tests pass, including concurrent separate database
connections, idempotency, invalid state transitions and unknown schema rejection.
Real child-process tests pass: independent owners, single supervisor lock,
capacity rejection, startup, forced child death, surviving sibling, explicit
stop and owned-directory cleanup. Repeated failed starts, startup timeout,
disconnected-run expiry, forced coordinator death/child death and durable
receipt-based recovery pass. Recovery never reconstructs a kill target from
a saved PID or claims to restore UI memory.

`serve_owned_ui` is the shared Bonito child adapter. Its bounded startup phase
prepares each registered renderer and closes probe sessions before reporting
ready. The initial browser fixture exposed cold workbench compilation delaying
sibling asset requests; moving UI-only renderer preparation into startup fixes
that boundary without increasing HTTP timeouts or importing scientific code.

The real-browser gate verifies two frames sharing one run, a second independent
run, binary Observable updates, disconnect/reconnect without DOM replacement,
four light/dark changes, existing X-ray metadata/CSS controls, surviving-host
interaction after sibling stop, shared-theme recovery HTML, and explicit restart
into a new empty run. These checks use isolated Chrome, not the user's session.

Final P1 logs: `runtime-unit-tests.log`, `p1-integration-final.log`, and
`bonito-gateway-final.log` in the baseline diagnostic directory. Runtime profile
records are passive; scientific executor hooks/preparation remain P4 work.

Commands so far:

```sh
julia --startup-file=no --threads=2 --project=playground/runtime playground/runtime/test/runtests.jl
julia --startup-file=no --threads=2 --project=playground/runtime playground/runtime/test/integration.jl
bash playground/runtime/test/run-bonito.sh
```

## Decisions and changes

- The approved plan is active. The opt-in `lcm runtime` launcher now publishes
  the shared catalogue and isolates explicitly started application UIs. The
  scientific showcase remains a passive skeleton until P6 installs its views.
- The initial ribbon failure read the prior hover frame immediately after CDP
  input acknowledgement. Added a rendering settle and hit-test diagnostics;
  a subsequent complete browser run passes without changing application CSS.
- New Julia API documentation follows the `julia-docstrings` skill, including
  field ownership/lifetime documentation and explicit error contracts.

## P2 implementation evidence

The single passive catalogue in `src/applications/Catalogue.jl` supplies both
the published JSON and the runtime registrations. It describes the scientific
showcase, template workbench, toolkit gallery, starter deck and hostile deck.
Catalogue category and entry surface are separate: a gallery can be a published
document whose original Bonito widget factories all belong to one isolated run.

The scientific entrance is `/`; full foundation references moved to `/dev/`.
Public deck/workbench selectors perform explicit, idempotent admission only.
Reading the catalogue, navigating static slides or browsing developer reference
pages allocates no processes. Uninstalled scientific adapters cannot launch.
The final scientific narrative and operations are still P6 work, not complete.

Published widget frames and live deck frames share their immutable owned-run
prefix. A standalone developer publisher retains its original widget routes;
the authenticated runtime never exposes the legacy unassigned job panel.
Original upload fields now derive their endpoint from their Bonito server and
carry the mutation marker through the gateway. No widget markup is replicated.
Child temporary uploads live beneath the owned run's scratch directory.

The end-to-end gate starts the actual `lcm runtime` CLI, checks 18 public/developer
routes in both themes and desktop/compact layouts, crawls generated navigation
links, launches a real deck/workbench/gallery, checks live-frame identity, then
uploads, downloads and removes a file through the owned namespace. It shuts the
CLI down with live WebSockets open and verifies terminal run records and empty
owned scratch. Log: `p2-catalogue-verified.log` (pass, including 3 teardown checks).

The expanded cleanup audit found and fixed a real CLI race: the main task could
exit after listener closure while background child cleanup was still yielding.
Shutdown is now joined through a shared lock. Teardown also reaches all owned
hosts even when one ownership receipt is damaged; unverified scratch is preserved
and reported for operator attention. Four earlier, finished test directories were
reconciled using their ownership receipts; no unrelated resources were removed.
Owned UI commands reuse existing package caches without launching detached
package-precompilation subprocesses.

Verification logs in the baseline directory:

- `p2-architecture-unit.log`: 132 runtime assertions, including dependency,
  credential, shared-style and passive-catalogue boundaries.
- `p2-final-runtime-integration.log`: 100 HTTP/WebSocket/process assertions,
  including damaged-receipt sibling cleanup and coordinator-death recovery.
- `p2-final-publisher.log`: all groups pass, including upload lifecycle and
  440 visual-contract assertions.
- `p2-final-presentation.log`: full presentation gate passes, including native
  fragments, equation notes, overview/static modes and linked-placeholder PDF.
- `p2-final-ribbon.log`: shared ribbon/toolbar, navigation, tab, segmented and
  data-row selection round trips pass across gallery/embedded/workbench hosts.

The new runtime environment is resolved separately; existing publisher, worker
and protocol manifests have not been replaced. No NATS broker or scientific
worker was required for the P2 browser gate.

## P3 work in progress

`protocol/src/Runtime.jl` adds passive profile advertisements, challenge-bound
worker reports, assignment fences, lease controls/acknowledgements and assigned
scientific requests. The nested scientific payload remains the validated v1
record; v2 ownership and targeted subjects are separate. Strict decoding rejects
unknown/duplicate nested fields, boolean-to-integer coercion, malformed identities,
unsupported versions, invalid lease durations and oversized messages.
`p3-protocol.log`: 43 new assertions plus all 43 existing v1 assertions pass.
Wire validation alone is **not** broker authentication or a granted lease.

`runtime/src/WorkerStore.jl` adds administrator-only pending enrollment from an
operator-provisioned identity binding, explicit approval/drain/disable, durable
idempotent responses and revision checks. Replaying an old approval returns its
original response without reversing a newer drain/disable. Registration is not
presence or readiness, and disabling does not implicitly kill existing work.

SQLite schema 2 adds workers and control-action bookkeeping. Existing schema 1
requires explicit `lcm runtime migrate`: the same supervisor lock excludes a live
coordinator, a private consistent snapshot precedes a transactional migration,
and repeated migration is inert. New databases initialize directly at schema 2.
TOML configuration remains version 1. No user database was migrated; tests use
temporary copies. `p3-enrollment-unit.log`: 164 assertions pass, including 32
enrollment/migration checks. `p3-registration-integration.log`: all 100 existing
gateway/process assertions still pass against the new schema.

The runtime-only manifest now references the existing protocol and vendored NATS
client. Broker credentials, live inventory, lease allocation, worker agents and
control-panel components are not yet connected. The next transport work must
enforce per-worker subjects and private reply inboxes at the broker, not trust a
payload worker ID or grant broad cross-worker JetStream access. SQLite remains
the allocation authority; no independent broker lease table is introduced.

Transport inspection found an unbounded initial NATS handshake. The maintained
fork now applies a finite deadline to TCP/INFO/TLS/authentication, closing only
the owned attempt socket. A strengthened fixture proves that the authentication
case actually sent CONNECT/PING before timing out (`p3-broker-deadline-confirmed.log`,
16 checks). The same inspection found that optional server TLS could bypass an
explicit client TLS requirement; encryption is now selected when either side
requires it. The added TLS-stall/no-plaintext fixture is being verified in
`p3-broker-handshake-security.log`. Real broker/TLS integration and the full
legacy worker regressions remain required after these transport changes.

### P3 verified checkpoint — control transport (not yet the complete control plane)

- Challenged in-memory presence is separate from durable registration. Reports
  require the current coordinator challenge, advancing sequence, a non-retired
  boot, approved profile fingerprints and non-inflated capacity. Public reports
  omit challenge/coordinator nonces. Restart begins with unknown presence.
- `BrokerPolicy.jl` defines fixed coordinator/worker identities, private reply
  prefixes and per-worker JetStream stream permissions. Worker credentials may
  pull only their own precreated consumer; they cannot create broader filters,
  impersonate another report subject, or issue lease grants.
- `BrokerTransport.jl` uses a dedicated control connection, per-worker bounded
  queues, fair nonblocking polling, strict passive decoding and no disconnected
  publication replay. Endpoint configuration is inert, rejects credential-bearing
  URLs and requires verified TLS except explicit literal-loopback development.
  Raw library diagnostics are excluded from the runtime's owned connection tasks.
- Operator stanzas contain environment-variable references, not secret values.
  Existing v1 broker configurations now scope consumer APIs/acknowledgements to
  `LCM_JOBS`; v1 credentials cannot use their former broad consumer API grant to
  access v2 streams. Existing streams/data are not deleted or remapped.

The real legacy recovery gate exposed two pre-existing live-path defects:
ambiguous `validate` after EzXML imports, and executor bootstrap time being charged
to a five-second calculation budget. Validation is now qualified. A finite
bootstrap acknowledgement precedes calculation timing; it explicitly does not
claim scientific preparation. The original calculation deadline is unchanged.

The maintained NATS fork additionally forwards actual ACK/NAK/TERM bodies,
cleans failed SUB/UNSUB/request resources, closes reply timers in finally, drops
uncertain buffered publications when replay is disabled, and removes credential-
bearing URLs from reconnect messages. These changes are enumerated in the fork's
`LCM_PATCHES.md`.

Verified logs under `/tmp/lcm-runtime-platform.9uEfiTvj/`:

- `p3-probe-protocol.log`: 44 v2 plus 43 existing v1 assertions.
- `p3-bootstrap-worker.log`: 78 worker assertions, including cancellation,
  independent startup and real LineCableModels scientific parity.
- `p3-bootstrap-nats-tls.log`: complete publisher suite, all 43 real broker
  lifecycle checks, 3 authorization and 3 mutual-TLS checks. Crash recovery
  passes without increasing the calculation deadline.
- `p3-control-unit.log`: 220 runtime assertions, including 24 inventory and
  32 broker-policy/endpoint checks.
- `p3-transport-cleanup.log`: 48 transport checks, including real stalled
  TCP/auth/TLS peers, no uncertain replay, and request/subscription cleanup.
- `p3-control-broker.log`: 26 actual TLS broker assertions with two independently
  credentialed worker connections, exact grant/report/ack targeting, spoof
  rejection, cross-worker and legacy-v1 pull denial, immediate NAK redelivery,
  TERM completion and teardown. Owned broker removed and temporary credentials
  deleted; non-secret broker diagnostics retained at
  `/tmp/lcm-runtime-control.QESm0o9c/`.

P3 remains incomplete: these tests exercise the transport with passive records,
not a usable lease manager or executing agent. Durable allocation, actual lease
lifetime enforcement, coordinator/agent loops, control-panel actions and reusable
diagnostic controls remain to be implemented. The full legacy TLS gate must be
repeated after the latest request-finally/no-replay changes; the dedicated v2
TLS gate already includes those changes.

### P3 durable allocation checkpoint

SQLite schema 3 adds immutable assignment fences and capacity bookkeeping, with
unique live run/role and owner/request constraints. Explicit backed-up migration
supports schema 1 and schema 2; approval and idempotent enrollment history survive.
No user database was migrated. Configuration TOML remains version 1.

`AssignmentManager` authorizes the actual run owner and declared role/profile,
requires current approved compatible presence, then admits in one IMMEDIATE
transaction. Automatic placement uses available capacity; pinned placement never
falls back; a dedicated reservation excludes other runs while allowing its own
sibling roles. Global, per-owner and per-run counts include pending and unreconciled
reservations. Operator allocations are still charged to the run's real owner.

These are **reservations**, not usable grants. Restarted inventory or a fresh
worker report alone does not erase old reservations or establish executor warmth.
Grant acknowledgement, lease expiry/release and actual agent supervision remain
the next implementation boundary.

- `p3-allocation-unit.log`: 263 runtime assertions pass, including 43 new
  ownership, placement, generation, capacity and schema-2 migration checks.
- `p3-control-runtime-integration.log`: all 148 gateway/process/transport checks
  pass against schema 3 (100 pre-existing integration checks plus 48 transport).

The small broker acceptance harness is explicitly unignored so the repeatable
test is included with the runtime sources.

### P3 lease-lifetime checkpoint

`AgentLeaseLedger` now enforces exact coordinator/worker boot, profile,
run/role/generation and revision fences. Its own monotonic clock controls
authority. Duplicate current commands return the same acknowledgement without
renewing expiry. Expiry, lost coordinator presence and release revoke new starts
immediately; closing slots remain occupied until the owning supervisor explicitly
confirms cleanup. Out-of-order release installs a generation tombstone, so a
delayed grant cannot recreate already released authority. Finite history fails
closed instead of forgetting fences. New coordinator probes retire old authority
before reconciliation, including delayed old probes.

`LeaseCoordinator` joins SQLite reservations to that protocol. Only an exact
acknowledged grant establishes temporary authority; persisted active rows do not.
Renewal preserves prior acknowledged authority while awaiting its bounded reply.
Late grant/renewal replies cannot revive revoked use. Missing replies and broker
outages initiate/retry the same release, holding capacity until cleanup is
confirmed. Abandoned unissued reservations have bounded grace, and stopped runs
initiate release. A freshly challenged replacement report may retire old
incarnations only under the agent's cleanup-before-announcement protocol.

Verified evidence:

- `p3-agent-lease-unit.log`: 318 runtime assertions, including 55 worker-side
  duplicate, expiry, cleanup, generation-history and coordinator-replacement checks.
- `p3-coordinated-lease-unit.log`: 373 runtime assertions, adding 55 joined
  acknowledgement, lost/reordered delivery, broker outage and restart checks.
- `p3-lease-lifetime-unit.log`: 382 assertions, adding abandoned-reservation and
  stopped-run release checks.
- `p3-lease-tls.log`: 42 real TLS checks: 26 identity/ACL/targeting checks and
  16 checks of the actual SQLite + coordinator + agent ledger grant/renew/release
  path. This fixture owns no executor; it explicitly does not claim physical
  process teardown has been implemented by the ledger.
  Broker removed and temporary credentials deleted; non-secret diagnostics:
  `/tmp/lcm-runtime-control.9W7pDIPY/`.
- `p3-final-legacy-tls.log`: after the final NATS fork request cleanup/no-replay
  changes, all publisher groups, 43 lifecycle checks, 3 broker authorization and
  3 mutual-TLS checks pass again. The prior recovery timeout is resolved.

**P3 is still in progress.** The coordinator/agent scheduling loops, authenticated
configuration/control-panel routes, reusable selectors/status/diagnostics,
targeted durable job/result orchestration and executing-agent lifecycle are not
yet wired. No new runtime control capability is advertised to browser consumers
until those boundaries are implemented and tested. Scientific profile preparation
and terminal isolation remain P4/P5, not claims made by lease acknowledgement.

Additional P3 evidence:

- `p3-allocation-process-race.log`: 387 runtime assertions, including a real
  two-Julia-process barrier race against one SQLite worker slot. Exactly one
  reservation succeeds; the other sees unavailable capacity.
- `p3-assigned-results-protocol.log`: 99 protocol assertions (56 v2 and the
  original 43 v1). `AssignedResult` wraps the existing scientific result with its
  assignment fence and rejects foreign worker/environment provenance, extra
  nested fields and coercible artifact/failure fields. Result subjects carry
  worker boot, lease, generation and job UUID.
- Worker broker permissions now permit lookup of only their own result stream,
  needed to recover result-before-ack without repeating an already completed
  operation. That small permission extension still needs the next real TLS
  job/result transport gate; earlier TLS results above predate that extension.

Next implementation work is the targeted durable job/result adapter, followed
by coordinator/agent loops and authenticated configuration/API/UI integration.
No result adapter or executing agent is implied by the new wire record.

### P3 targeted job/result transport checkpoint

`AssignedJobs.jl` adds a separate scientific transport connection. Provisioning
creates bounded per-worker v2 job/result streams and one fixed worker consumer.
Existing configuration must match; it is not silently resized, purged or replaced.
Jobs have work-queue retention, finite message/byte/age limits and three maximum
deliveries. Results have bounded retention and one immutable outcome per subject.
The legacy v1 streams remain separate and unchanged.

Publication requires the owned acknowledged fence and an operation permitted by
that profile. The worker performs one finite pull and rejects malformed, foreign,
expired or incompatible authority before returning any request to an executor.
Result persistence matches the exact request, stores before input ACK, and returns
the already-stored original on retry instead of exposing a recomputed replacement.
Job/result deduplication IDs are scoped to the assignment. This adapter performs
no scientific execution or implicit preparation; the fixed execution orchestrator
must still own readiness, queue admission, execution and terminal outcomes.

- `p3-job-stream-provision.log`: 44 real TLS checks, including repeated inert
  provisioning and the new worker-own-result permissions.
- `p3-assigned-jobs-tls.log`: 75 real TLS checks: 28 identity/provisioning, 16
  lease round-trip and 31 targeted durable job/result checks. Synthetic scientific
  results explicitly stand in for an executor in this transport-only fixture.
  Verified duplicate publication, own-stream lookup, cross-worker denial, stored
  result with lost input ACK, original-result preservation, malformed job
  termination and expired-lease rejection. Owned broker removed; test credentials
  deleted; non-secret diagnostics: `/tmp/lcm-runtime-control.XHjXv3kZ/`.
- `p3-assigned-runtime-unit.log`: all 387 runtime assertions still pass,
  including the independent-process SQLite allocation race.
- `git diff --check` remains clean.

P3 still needs the coordinator/agent service loops and configured/protected API
and component integration. The scientific execution orchestrator, preparation and
physical executor cleanup remain P4 work. No browser capability has been enabled
on the strength of synthetic-result tests alone.

### P3 configured services, protected API and shared client checkpoint

The optional runtime `[control] config_file` now loads strict server-owned
broker/profile/trust configuration without connecting or importing engines.
`ControlService` schedules finite background connections, challenged reports,
lease renewal and expiry independently of the public gateway. Existing persisted
trust cannot be silently expanded by changing configuration. Protected HTTP
endpoints share the same principal/origin/CSRF checks as owned UI runs; enrollment
looks up a provisioned ID rather than accepting browser credentials or launch
arguments. Registration revision/idempotency, role placement and explicit release
are enforced through the existing root methods. Response writes do not hold the
coordinator lock.

Structured diagnostics are bounded and owner-filtered, accept only fixed codes
and validated identities, and carry a fresh buffer epoch as well as a sequence.
This detects coordinator restart even when a replacement buffer's sequence has
already surpassed the browser's old cursor. Duplicate JSON fields and excessive
nesting are rejected at the shared gateway reader.

`AgentService` adds the matching worker scheduler and a checked
`AbstractAgentResources` contract. Recovery precedes the first announcement;
cleanup runs separately from heartbeats. Release/replacement acknowledgement is
withheld until the exact resource owner confirms completion. Failed cleanup
retains occupied capacity and retries at a bounded rate. Concurrent close calls
join complete teardown; unresolved physical cleanup remains an explicit error,
including a repeated close. An agent may advertise zero installed profiles so a
healthy supervisor remains diagnosable without becoming eligible for execution.
This is a compatible v2 payload relaxation, not a new execution permission.

Added `lcm runtime permissions` (read-only ACL/environment references),
`provision` (explicit bounded v2 stream creation/verification), and `check-agent`
(passive agent configuration validation). Documented their exact supported schema
and API in `runtime/CONFIGURATION.md`. Native/container launch still needs the
P4 concrete resource owner; there is intentionally no pretend runnable agent CLI
with an empty executor behind it.

Verified evidence:

- `p3-service-final-unit.log`: **486** runtime assertions, including strict
  control/agent config, CLI redaction, event epochs, required cleanup hooks,
  nonblocking cleanup, bounded retry and joined teardown.
- `p3-agent-service-tls-verified.log`: **119** assertions over the isolated
  real TLS broker: 28 ACL/provisioning, 16 lease round-trip, 31 targeted durable
  jobs/results, 20 coordinator scheduling and 24 joined production
  coordinator/agent recovery checks. The resource fixture explicitly owns no
  executor; this is not a scientific execution or container-isolation claim.
  Final event-epoch/client-asset changes are covered by the later unit/HTTP gates.
  Non-secret broker diagnostics: `/tmp/lcm-runtime-control.H6I80CVJ/`.
- `p3-control-runtime-full.log`: **219** integration assertions, retaining
  existing UI-host, ownership, shutdown and transport checks plus the 55 control
  HTTP and 16 stalled-broker checks.
- `p3-service-final-http.log`: those **71** control/stalled-broker HTTP checks
  pass after the final service/event changes. Standalone cold HTTP specialization
  took approximately 5.5 seconds (`p3-control-offline-recheck.log`); repeated
  health calls while INFO was deliberately withheld stayed below one second.
  The combined already-partly-compiled fixture measured 1.8 seconds for its
  first health specialization. Cold compilation is not reported as steady-state
  responsiveness.
- All isolated `lcm.runtime.test=control-v2` containers were removed, credentials
  deleted by the fixture, and `git diff --check` is clean. Pre-existing containers
  remain untouched.

The shared browser client now lives in `assets/runtime-client.js`, exposed as
the passive `/runtime/assets/runtime-client.js` resource. Import/construction is
inert. It coalesces polling per run/window, marks retained data stale on failure,
separates log polling from inventory, bounds responses/history, retains event
epochs, and never automatically replays an uncertain mutation. Post-action
reconciliation cannot reuse an older in-flight snapshot. Last-subscriber teardown
aborts requests and removes timers/listeners.

`p3-client-asset-http.log` passes **57** protected-control/asset assertions,
including the served shared client and unchanged authorization boundaries.

`node test/integration/runtime_client.mjs` passes the client transport/lifecycle
fixtures: passive construction, coalescing, stale data, independent stalled logs,
bounded history, restart gaps, visibility, no mutation replay, exact explicit
retry identity, fresh reconciliation and teardown. This is a client unit gate,
not a browser visual-conformance result. The actual WorkerSelector,
PreparationStatus, WorkerDiagnostics and control-panel composition still need
their common render paths, CSS/X-ray ownership and gallery/workbench/deck tests.

**P3 remains in progress.** Remaining P3 work includes those components/panel,
remote execution-stage/event correlation and runnable configured agent integration
with the P4 resource supervisor. P4 preparation/physical cleanup, P5 private
terminal, P6 scientific migration and P7 aggregate release gates remain open.

### P3 continuation — shared runtime controls and protected panel

Implemented `RuntimeClient`, `WorkerSelector`, `PreparationStatus`,
`WorkerDiagnostics`, and `WorkerControlPanel` in
`src/widgets/RuntimeControls.jl`. The Julia values are passive run/role metadata;
the browser owns subscriptions and the coordinator remains the authorization
and assignment authority. They compose as ordinary Bonito DOM values and expose
source-scoped X-ray metadata without credentials, event payloads or Julia code
evaluation. Their new API docstrings follow the Julia docstring skill.

One `assets/runtime-controls.js` renderer is used by Bonito wrappers and the
protected `/runtime/control` page. Its stylesheet owns layout only and consumes
the established brand/native-form tokens. The panel's optional `?run=UUID`
selection is authorized before rendering, and its selectors come from that
application's registered requirements. Merely opening it allocates nothing.
The administrator may enroll only provisioned identities and apply explicit
revision-checked approval/draining/disabled changes. The run list is an explicitly
labelled page snapshot; worker status and bounded events update live.

Selectors preserve unsent choices and keyboard focus through polling. An offline,
removed, full or incompatible pinned worker stays visible but cannot be submitted;
the UI never chooses a fallback. An uncertain mutation is not replayed
automatically: its explicit retry reuses the exact original request/input, while
other control mutations remain disabled until that uncertainty is addressed.
Assignment/release, registration and preparation remain distinct. Preparation
still reports unknown because P4 has not supplied executor preparation evidence.
Diagnostics currently show the bounded registered control events, not arbitrary
worker logs or terminal output.

The shared client now handles a valid owned UI run whose publisher has worker
control disabled, instead of turning the unavailable assignment endpoint into a
misleading connection failure. Removing the final component closes its client,
requests, timers and observer subscription. The gallery documents and mounts all
four views (with inventory-only selector specimens and no scientific allocation).
Developer navigation and run recovery pages link to the protected control panel.
`runtime/CONFIGURATION.md` documents the actual APIs and their remaining limits.

Verified evidence (all commands exited successfully):

- `p3-controls-full-publisher.log`: **1,057** assertions across all 23 publisher
  test sets, including **33** new runtime-component and **447** visual-contract
  assertions. Route/source audit includes the new gallery entry.
- `p3-controls-full-runtime.log`: **486** assertions across all 33 lightweight
  runtime test sets; no scientific or Bonito import was added to that environment.
- `p3-controls-final-http.log`: **68** actual HTTP assertions for control APIs,
  assets and protected panel selection, including anonymous/foreign-run rejection,
  duplicate query rejection, role derivation and zero UI-host allocation on render.
- `node test/integration/runtime_client.mjs`: shared client transport/lifecycle
  fixtures pass, including the disabled-control/owned-run case.
- `p3-controls-final-browser.log`: real Chrome exercises the actual shared client
  and renderer against an explicitly labelled local HTTP state fixture. Verified
  draft/focus retention, offline pinned state, light/dark/light transitions,
  native option contrast, 1280/390-pixel layouts, exact explicit retry, release,
  revision-bearing registration, independent event failure, stale data and final
  subscription teardown. Screenshots and isolated browser diagnostics:
  `/tmp/lcm-runtime-controls.4rUp3E/`. This fixture is not a broker/executor test.
- `p3-controls-bonito.log`: actual gateway + separate Bonito UI processes mount
  the same controls in gallery and workbench frames, preserve theme parity and
  X-ray source inspection, and retain the existing binary sockets, reconnect,
  cross-run isolation, owned stop, surviving peer and clean restart gates.
  Diagnostics: `/tmp/lcm-runtime-bonito.ZqAnng4X/`. Its owned host scratch is empty
  after joined teardown and no fixture processes remain. The fixture now uses
  an explicit stdin shutdown channel and finite shell cleanup instead of an
  unbounded signal/wait sequence.
- `p3-controls-build.log`: the complete Quarto publisher rebuild succeeded;
  generated developer/gallery links and the new frame route were checked.
- `git diff --check` is clean. No containers were started for this checkpoint;
  pre-existing containers and unrelated user data remain untouched.

The first browser fixture attempt exposed a redeclared test variable; the
corrected final browser run passes. Integration review also identified the
disabled-control/owned-run presentation issue above, which now has a client
regression test. These fixes do not constitute scientific readiness.

**P3 remains in progress.** The shared controls/panel implementation is now in
place. Next work is remote executor preparation/stage/event correlation and
the P4 concrete native/container resource owner plus runnable agent CLI.
P5 private terminals, P6 scientific showcase/workbench integration, and P7 full
release/recovery/deployment verification remain open. Physical Docker and CPU
cgroup verification remain the previously recorded external constraints; they
have not been bypassed or reported as successful.

### P4 continuation — separate scientific profiles and preparation core

Added the lightweight `worker/core/` package and pinned, independent
`worker/profiles/line-parameters/` and `worker/profiles/power-flow/` environments.
The core imports neither numerical engine, broker, artifact client nor UI package.
The line environment contains LineCableModels but not PowerImpedance; the power
environment contains PowerImpedance but not LineCableModels. Neither contains
Bonito/NATS/AWS/AWSS3. Their manifests were resolved offline from installed packages.
The PowerImpedance source revision remains the existing approved revision.

The core and legacy worker reuse the same existing `OperationRegistry.jl`,
`Cache.jl`, `Executor.jl`, and scientific operation-adapter sources. No scientific
operation implementation was copied into a parallel path. The new core package
requires the monorepo worker source tree when packaged, as documented in
`worker/profiles/README.md`. The legacy v1 combined registry is preserved.

Four checked `AbstractScientificProfile` hooks own registry, preparation input
validation, representative work and cleanup. `PreparedWorkload` records evidence
and dependent model keys. Cached evidence is invalidated if those models expire
or are evicted. Line preparation runs the actual geometry at endpoint frequencies;
power-flow preparation actually constructs, solves and linearizes the existing
OHL/UGC model. This is representative preparation, not a promise that all future
Julia specializations are compiled.

Extended the existing `ExecutorSupervisor` with operator-owned profile launch
commands and explicit `prepare_supervised!` on the same framed channel. Bootstrap
remains distinct from preparation. EOF and oversized input/output lines are
bounded, operation elapsed-time checks use the monotonic clock, and a failed
progress/log observer retires the process so leftover frames cannot become the
next operation's result. Oversized commands are rejected before process creation.
The existing cancellation/deadline path remains shared.

Preparation caches now enforce admission while builders are pending, bounded
retained bytes, monotonic idle expiry, and cleanup identity invalidation. A late
builder cannot refill a cleared cache or overwrite a newer same-key entry.
The byte limit is estimated retained-value accounting, not a claim about peak
solver/package/process memory. Cached failures rethrow their stored cause instead
of manufacturing a new builder; the legacy test now documents that distinction
from a newly failed TaskFailedException. API docstrings were added with the Julia
docstring skill and the actual ownership/resource limits are documented.

Verified evidence:

- `p4-core-import.log`: successful core import with assertions excluding both
  numerical engines, NATS and Bonito. The initial attempt exposed that the installed
  RequiredInterfaces macro also needs an alias for a parameterized argument type;
  the corrected required declarations and negative missing-hook tests pass.
- `p4-core-final-unit.log`: **39** assertions covering import boundary, bounded
  framing, rejection before launch, cache admission/bytes/expiry, single flight,
  cleanup invalidation, required hooks, and evidence whose model was evicted.
- `p4-legacy-final-worker.log`: **80** assertions retaining legacy operation,
  cache/artifact, cancellation, engine-isolation and LineCableModels numerical
  parity gates, plus the failed-observer retirement regression.
- `p4-two-profile-processes.log`: **24** assertions over actual separate native
  Julia children. Real line preparation reports miss then hit, performs a real
  frequency scan, rejects a power-flow operation, and returns to miss after
  process replacement. Power-flow preparation solves and linearizes the real
  case; an independent line request completes while that preparation runs, and
  stopping power flow leaves line execution usable. Parent import boundaries
  and profile manifests are checked. The final observer-failure/oversized-command
  guards added after this numerical run are covered by the later core and legacy
  gates above; P7 retains the full aggregate rerun.
- Measured line representative workload: **7.94 s** cold and **0.027 s** cached
  (not package startup). Independent line request: **0.130 s**. Repeated power-flow
  preparation: **0.027 s**. The entire two-profile test set, including startup and
  preparation, took approximately **3m17s**. These are local measurements, not
  latency guarantees or a complete power-flow cold-start breakdown.
- All launched numerical test processes exited, `git diff --check` is clean,
  and no containers or unrelated services were touched. Existing duplicate engine
  docstring/precompilation warnings remain visible; no engine source was altered
  to suppress them.

**P3/P4 remain in progress.** These local process tests do not yet implement the
agent's lease-owned resource manager, full source-fingerprint preflight, durable
preparation/job request orchestration, executor-correlated remote stage events,
native/container resource policy, or runnable agent CLI. The gateway still
advertises `assigned_execution=false`; reusable preparation controls still report
unknown. Next work connects the tested core into that authoritative agent/lease
path, without importing engines into the runtime coordinator. P5–P7 remain open
and the original external CPU-controller/Docker constraints are unchanged.

### P4 continuation — source identity and bounded process retirement

Added `runtime/src/EnvironmentFingerprint.jl` and the read-only
`lcm runtime fingerprint --project DIRECTORY` action. Both scientific packages
can be inspected without engine imports, child processes, network access or
broker configuration. `verify_native_environment(profile)` recomputes and checks
the approved digest; it does not advertise readiness or perform isolation preflight.

The bounded digest covers the selected version-specific/general manifest,
root package and every mutable manifest path dependency's Project.toml, src/,
ext/, deps/, artifact declarations and local preferences. Three small
`RuntimeSources.toml` declarations include the shared executor/operation files
which intentionally live outside the package source directories. Tests derive
literal external includes and require corresponding declarations. Nonliteral
includes/external data still require explicit author declarations; this is not
a Julia static-analysis framework. The reader rejects source/UUID/version drift,
unidentified unpinned dependencies, spoofed stdlib identity, alternate unsupported
project filenames, symlinks, excessive entries/bytes/depth, and concurrent source
changes observed over the complete inspection. Third-party Git tree pins identify
dependencies but do not attest an operator's mutable installed depot or artifacts.

The same existing scientific `ExecutorSupervisor` now owns bounded IO tasks.
Submitting to a child that never consumes stdin remains deadline/cancellation
aware. Retirement attempts the shutdown frame asynchronously, escalates through
finite TERM/KILL waits, closes pipes and joins owned tasks. A failed physical/pipe
cleanup retains its original process handle and cannot falsely release a slot.
Repeated stop joins the same cleanup; a malformed response missing its type now
retires the process too. Already-canceled startup is rejected before allocation.
The operation lock and independent stop lock preserve the ability to retire an
occupied executor without waiting for its numerical response.

The primitive still owns one trusted native process, not a subprocess tree or
container sandbox. Its future agent adapter must enforce those stronger resource
boundaries and call it under exact lease authority. Julia API documentation and
operator documentation were written using the Julia-docstring skill to state
these limits explicitly, rather than label source verification as preparation.

Verified evidence under `/tmp/lcm-runtime-platform.9uEfiTvj/`:

- `p4-fingerprint-final-runtime.log`: **547** assertions (the previous 486 plus
  61 source-inspection/CLI/conformance assertions), all passing. The intermediate
  fingerprint runtime gate passed 540 before the additional identity/include tests.
- `p4-fingerprint-cli.log`: the actual `./lcm` launcher prints one source digest
  for the line profile. Source changes correctly invalidate that digest; do not
  copy a historical test digest into a deployment configuration.
- `p4-bounded-executor-final-pass.log`: **75** core assertions, including real
  full-pipe submission timeout/cancel, scheduler progress, confirmed SIGSTOP then
  KILL escalation, malformed response, retained unresolved pipe cleanup and joined
  repeated stop, plus the prior framing/cache/preparation evidence gates.
- The initial libc SIG_IGN fault fixture was unreliable with Julia's signal-wait
  thread. The final test explicitly suspends its original Process handle, verifies
  Linux reports it stopped, then requires KILL retirement. Julia does not expose
  SIGSTOP and this host's external `kill -l` only converts numbers to names; the
  final fixture obtains the number from the verified Bash builtin. Failed fixture
  attempts remain in diagnostic logs. No shutdown requirement was relaxed. Expected
  Julia TERM diagnostic stacks from deliberate non-cooperating children remain
  visible in the test log, not misclassified as unexpected suite failures.
- `p4-bounded-executor-legacy.log`: **80** assertions, including real
  LineCableModels numerical parity, existing v1 supervision/cancellation and caches.
- `p4-fenced-core-two-profiles.log`: **24** real separate-profile assertions pass
  after the shared supervisor changes. Despite this diagnostic filename, this is
  a local core/process test, **not** an end-to-end lease-fencing integration gate.
  Representative line workload: 8.16 s cold, 0.053 s cached; independent line
  request during power preparation: 0.131 s; repeated solved/linearized power
  preparation: 0.052 s. The power-independence test set including its process
  startups completed in approximately 3m09s. These are measurements, not guarantees.
- All owned numerical and fault-injection children exited. `git diff --check`
  is clean. No containers, pre-existing services or unrelated user files changed.

**P3/P4 remain in progress.** Source verification is implemented but not yet
connected to a concrete agent resource owner/launch preflight. That owner, exact
lease/executor-correlated preparation and job orchestration, physical resource
policy/container adapters and runnable agent CLI remain the next work. The gateway
still correctly advertises `assigned_execution=false`; heartbeat and assignment
controls do not imply a prepared scientific process. P5–P7 and the previously
recorded actual-Docker/CPU-controller verification constraints remain unchanged.

## P4 — lease-owned scientific execution checkpoint

`ScientificResources` now binds once to the agent's exact live lease ledger and
owns one independently supervised process/cache per assignment. Required physical
driver hooks cannot replace its admission, source/lease checks, cancellation or
result checks. Registration remains passive; only explicit preparation admits a
new process. Concurrent identical preparation shares one task, different work is
rejected while busy, and the most recent identical completed job retry reuses its
task. Durable replay/result-before-ack remains the broker's responsibility.

Preparation is measured in the actual child: the shared execution core now
inspects retained model dependencies without rebuilding or extending their TTL.
An assigned operation requires that evidence and cannot silently start a cold
replacement. Activity and preparation are separate status dimensions. Only
bounded progress/counts and fixed failure codes cross this owner; raw engine
messages and private exception contexts do not. Loss of lease authority cancels
work independently of renewal/control scheduling. Unresolved physical cleanup
retains capacity; an already-closed agent can retry that same teardown without
reopening admission.

Verified logs under `/tmp/lcm-runtime-platform.9uEfiTvj/`:

- `p4-owned-recovery-final-runtime.log`: **617** assertions pass, including
  independent preparation, lease expiry, exact-request cancellation, cache eviction,
  source/preflight rejection, private error context and eventual cleanup retry.
- `p4-owned-core.log`: **86** assertions pass, including child inspection that
  neither extends TTL nor launches a missing process, plus bounded IO retirement.
- `p4-owned-legacy-worker.log`: **80** assertions pass, including numerical parity.
- `p4-owned-real-profiles.log`: **21** assertions pass with real line-parameter
  and power-flow engines in separate processes under live monotonic leases and
  independent renewal tasks. A line request completed during power preparation
  (0.732 s); repeated power preparation reused its process/models (0.137 s).
  Releasing power flow left line execution usable. These are measurements, not
  timing guarantees. All owned children exited.

The real-engine gate uses an explicit **test-only native driver**, not a
production OS-resource preflight or a broker job-delivery integration. Production
native/container drivers, runnable agent CLI, remote preparation/job orchestration,
executor-stage transport and protected API/client integration still remain. The
gateway continues to advertise `assigned_execution=false`. No containers were
started or existing services modified. P3/P4 stay in progress; P5–P7 retain their
full scope and the external CPU-controller/actual-Docker constraints remain.

## P4 — shared host-command and container prerequisite checkpoint

Container-engine discovery now lives in one dependency-free shared source,
`common/container_engine.jl`, consumed by the existing Compose CLI and the
lightweight runtime. Existing Docker-first preference, explicit selection and
Podman-shim recognition are retained; agent discovery does not require Compose.
The runtime declares this external source in its fingerprint inventory.

`CommandRunner` supplies the physical adapter's bounded host-command ownership:
finite admission, combined stdout/stderr byte limits, independent cancellation,
deadlines, TERM/KILL retirement and joined readers. A failed physical/pipe cleanup
retains its original handle until retry succeeds. Private output and original
exception contexts do not enter the fixed error projection. The container CLI
receives only operator configuration references, not broker/proxy/storage secrets
or inherited remote-engine selection variables. It does not receive user input.

`lcm runtime check-host [--runtime auto|podman|docker]` is now a usable read-only
operator diagnostic. It rejects missing required engine prerequisites, remote
engine contexts and unverified schemas. Docker's local socket is inspected and
returned as an explicit command prefix; Podman is explicitly local. This is
engine prerequisite evidence, **not** per-container effective-limit attestation.
No profile is advertised or made ready by this command.

Verified evidence under `/tmp/lcm-runtime-platform.9uEfiTvj/`:

- `p4-container-command-tests.log`: **100** focused assertions pass for command
  timeout, output flood, exact cancellation, independent work, full admission,
  unresolved cleanup/retry, private environment filtering, and both engine schema
  fixtures. The initial smoke test exposed an unsupported `readbytes!` keyword
  on Julia Pipe; it was corrected before the passing focused and full gates.
- `p4-container-foundation-runtime.log`: **717** assertions pass, exit 0,
  including all lease-owned preparation/cleanup regressions.
- `p4-container-foundation-publisher.log`: **1,057** assertions pass, exit 0,
  including unchanged Compose resolution/construction, registered components,
  visual contracts, workbench and X-ray contracts.
- `p4-container-host-preflight.log`: the actual `./lcm` command selects native
  rootless Podman and exits **2** with `cpu_controller_missing`, as required.
  Actual Podman info confirms only memory and pids are delegated. No actual
  Docker Engine is present. No containers were created or modified.
- All owned fixture/test children exited; `git diff --check` is clean.

**P3/P4 remain in progress.** This checkpoint supplies shared discovery, bounded
host commands and actionable prerequisites for the physical drivers; it does not
substitute those prerequisites for launch/receipt/kernel-limit verification.
Production native/container launch and recovery, runnable agent service wiring,
remote preparation/job orchestration and protected API/client integration are
still required. P5–P7 and the full approved finish line remain unchanged.

## P4 — exact physical ownership and container recovery checkpoint

`ResourceJournal` now records an acquisition intent before physical creation and
atomically binds its full runtime identity afterwards. A private directory,
single-writer kernel lock, bounded strict JSON, file/directory synchronization and
atomic replacement protect those records. Changed directories/locks, linked or
non-regular files, foreign ownership, duplicate leases and malformed data prevent
recovery actions. Opening is nonblocking even for a substituted FIFO. Validated
startup can remove incomplete-write fragments, but not unknown or corrupt data.
A killed agent releases its kernel lock while leaving recoverable receipts.

Receipts bind the full assignment fence, generated UUID name, backend and
**supervisor scope**. Scope identifies the inspected local engine, not merely a
`docker` or `podman` command. Docker daemon replacement, Podman graph-root
replacement and changed contexts cannot establish absence of old resources.
The Julia-docstring skill was used to document the distinction between cleanup
identity, live lease authority and actual preparation evidence at the API itself.

The actual `remove_owned_container!`/`recover_containers!` adapters use the
shared bounded command owner. An unbound intent resolves through exact name and
ownership labels, then commits the full ID. Scope is rechecked before mutations;
only that full ID can be stopped or removed. Failed inspection, failed removal,
malformed inventory or a changed engine retains unresolved ownership. A working
engine must confirm absence before forgetting a receipt. Already-removed targets
can be retried, but an externally restored, unrecorded container cannot be deleted
by replaying an old receipt. Independent receipts are still attempted when one
cannot be reconciled. Other backends remain for their own required adapters.

Verified evidence under `/tmp/lcm-runtime-platform.9uEfiTvj/`:

- `p4-container-recovery-tests.log`: **154** focused assertions pass before
  the additional three FIFO checks, covering journal ownership/death/recovery,
  exact scope/ID binding and both engine command fixtures.
- `p4-container-recovery-runtime.log`: **874** assertions pass, exit 0,
  including all **157** new journal/recovery assertions and prior control,
  preparation, authority, identity and cleanup regressions.
- `p4-container-scope-local.log`: **8** assertions pass against actual native
  Podman. Its inspected storage scope is stable, empty-journal recovery allocates
  nothing, and the complete container inventory is unchanged.
- `p4-container-stopped-recovery-final.log`: **13** assertions pass, exit 0,
  against actual Podman using the already cached digest-pinned official Julia
  image. One newly created stopped container retained every configured CPU,
  memory and PID limit, never ran a process (`Running=false`, `Pid=0`), and
  was recovered from its unbound receipt and removed by exact ID. The original
  inventory is unchanged and no `lcm-exec-*` test containers remain.
- The first stopped-container trial completed those 13 assertions and removed
  its own container, but its test wrapper had a Julia global soft-scope teardown
  error. The wrapper was moved into a function and the complete gate rerun
  successfully. A final audit also caught Julia's default exit cleanup of the
  diagnostic temporary directory; the test now explicitly uses `cleanup=false`
  so an unresolved future resource journal cannot be discarded at process exit.
  The initial trial had already removed its only resource before that directory
  disappeared. No existing user container was stopped or removed.
- `p4-container-recovery-preserve-failure.log`: an explicit missing-local-image
  negative gate exits 1 before acquisition, performs no pull/create, and leaves
  its private mode-0700 directory present after process exit. This confirms the
  corrected failure-artifact retention behavior; it is not an unexplained suite
  failure or an outstanding container.
- `git diff --check` is clean. All owned test children exited. The real
  container tests created only stopped acquisition metadata and removed those
  newly owned containers; **no container process was started**.

**P3/P4 remain in progress.** Receipt/container recovery is now implemented and
has real stopped-resource evidence. This is not a running-container quota,
scientific launch or terminal-isolation gate. Effective kernel enforcement,
production native/container launch, native-unit recovery, managed agent startup,
remote preparation/job orchestration and protected API/client integration remain.
Actual Docker execution is still unverified, this host still lacks delegated CPU,
and P5–P7 retain their complete approved scope.

## P4 — shared container policy, effective entry checks and image layout

`ContainerPolicy` now derives both scientific and terminal configuration from the
same approved profile and ownership receipt. `create_owned_container!` checks
host prerequisites, exact engine scope, installed image digest/contract, fixed
entry command and environment, and stopped-container policy. It creates neither
an implicit image pull nor a running process. Failed acquisition keeps its
receipt until physical recovery proves absence. Docker/Podman differences are
normalized inside this adapter, not in individual profile implementations.

The lightweight core now owns `ContainerLimits` and
`verify_container_isolation`. The fixed image entry guard checks effective Linux
UID/GID/groups, capability sets, seccomp, no-new-privileges, private cgroups,
finite CPU/memory/task limits, zero additional swap, network interfaces, approved
mounts, bounded tmpfs capacity and descriptor/message-queue/core-file limits.
It exits with fixed diagnostics before numerical imports or user evaluation on
missing or incompatible evidence. Configured limits, effective kernel evidence,
live lease authority and preparation are explicitly separate API contracts; the
Julia-docstring skill was used to retain those distinctions in the implementation.

There are separate immutable image targets for line parameters, power flow and a
lightweight Julia terminal. The terminal project has an offline-resolved pinned
manifest and no numerical/UI/broker dependency. Recipes copy only declared
source/lockfiles. The active project is a real directory at the original profile
depth: a local Julia probe exposed broken dependency/include resolution through
the initial short symlink, so that layout was replaced before acceptance. A
conformance test reconstructs the actual COPY instructions, verifies relative
dependencies and compares every image command with the shared launch policy.
Image builds and image provisioning remain explicit operator actions.

Verified evidence under `/tmp/lcm-runtime-platform.9uEfiTvj/`:

- `p4-container-policy-runtime-final.log`: **1,124** assertions pass, exit 0.
  This includes the prior 874 plus **175** shared-policy/acquisition assertions
  and **75** image-layout/real relocated Julia-import assertions.
- `p4-container-policy-core-final.log`: **149** assertions pass, exit 0,
  including **49** effective-kernel fixture assertions and **14** actual
  child-entry denial assertions. The evaluation marker is never executed after
  missing/host isolation; diagnostics are fixed and do not expose stack traces.
- `p4-container-stopped-policy-current.log`: **12** assertions pass, exit 0,
  against actual Podman with the final commands. One new stopped scientific
  resource and one new stopped TTY resource were inspected and removed by exact
  receipt-verified IDs. Both stayed `Running=false`, `Pid=0`. No process started.
  Earlier trials exposed Podman's cgroup/capability/rlimit/tmpfs inspection
  differences; those trials also removed their own resources. Their private
  inspection diagnostics remain in `/tmp/lcm-stopped-policy-WqUs8P` and
  `/tmp/lcm-stopped-policy-16jAjF`; neither owns an outstanding container.
- `p4-container-image-layout.log`: the standalone **75**-assertion layout gate
  passes with actual core imports and no scientific/UI/broker modules loaded.
  This is source/command conformance, not an OCI image build.
- `p4-terminal-environment.log`: terminal dependency resolution completed with
  `JULIA_PKG_OFFLINE=true`; no dependency download or image pull occurred.
- `p4-container-policy-host-final.log`: the real CLI still exits **2** with
  `cpu_controller_missing`. The rootless host does not have delegated CPU.
- All owned test child processes have exited, all newly created test containers
  were removed, and the original 16-container inventory remains unchanged.
  No `lcm-exec-*` container remains. `git diff --check` is clean. SIGTERM traces
  in the complete suites come from their deliberately terminated child fixtures,
  not a failed parent suite.

**P3/P4 remain in progress.** Checked stopped acquisition, effective entry checks
and image recipes are implemented; no running-container quota or prepared-model
claim follows from them. Next are managed physical scientific launch, native-unit
ownership/recovery and agent-death teardown, followed by runnable agent startup
and remote preparation/job API/client wiring. Actual image builds, effective
running-container checks and Docker integration remain unverified. P5 terminal
transport/UI, P6 scientific consumers and P7 release acceptance retain their full
approved scope. The goal is active, not complete, and mandatory limits have not
been weakened to fit this host.

## P4 — managed agent and scientific container-driver checkpoint

`ManagedAgentIdentity` and `agent_service_unit` now define the production
user-service boundary. `verify_managed_agent` reads typed systemd properties and
the actual process cgroup: exact worker unit, main PID, invocation, fixed Julia
start/recovery commands, complete-cgroup termination, no remaining-active state
after exit, and finite stop policy. A foreground caller or environment flag is
not accepted. The service rejects overrides that would retain a dead agent or
change post-stop cleanup semantics. The renderer is inert; installation and host
controller delegation remain operator actions.

The existing `lcm runtime` router now exposes `agent-unit`, `start-agent`, and
`recover-agent`, without another top-level CLI. Recovery captures the original
absolute journal and worker when rendering the unit and does not reread a mutable
agent config or broker credential after process death. A live kernel journal
lock rejects competing cleanup. Successful recovery restores neither lease
authority nor preparation; unresolved receipts remain owned. Native-unit receipts
still fail closed until that adapter exists.

`ContainerScientificDriver` implements the six required physical scientific
hooks. It admits only the installed/verified subset of approved scientific
container definitions and reports unsupported kinds/limits as unavailable. No
container is pulled or prepared during startup. Before physical acquisition and
execution it rechecks the service incarnation, engine scope, image and shared
container policy. Partial acquisition retains its receipt and slot. The attached
CLI remains one original supervised process; a narrow `discard_stderr` option
suppresses private engine diagnostics without permitting arbitrary subprocess
pipelines. Physical removal is attempted even if CLI retirement fails; ownership
is released only after both physical cleanup and joined process/pipe cleanup.

`serve_agent` composes this driver, `ScientificResources`, and `AgentService`.
Recovery precedes announcement; control/expiry remain independent of numerical
processes. This host starts a control-only agent with zero eligible profiles and
the honest `cpu_controller_missing` reason. Native and terminal profiles remain
explicitly unavailable in this driver. The Julia-docstring skill was used to
document the distinction between service lifetime, physical isolation, lease
authority and actual model readiness.

Final evidence under `/tmp/lcm-runtime-platform.9uEfiTvj/`:

- `p4-managed-runtime-narrow.log`: **1,202** assertions pass, exit 0,
  including **53** managed-service/config/recovery assertions and **25**
  concrete driver rejection/partial-ownership assertions.
- `p4-managed-core-narrow.log`: **155** assertions pass, exit 0, including
  **6** checks of the actual single-process redirected command, pipe cleanup
  and rejection of arbitrary command pipelines.
- `p4-managed-legacy-worker-narrow.log`: **80** assertions pass, exit 0,
  including real supervised engine/scientific parity regressions. Deliberate
  child-termination traces in these suites are not parent test failures.
- `p4-managed-agent-systemd-final.log`: **19** real lifecycle assertions pass
  (plus the separately counted 53 unit assertions), exit 0, on actual user
  systemd 252. The test rejected foreground startup, seeded one owned stopped
  Podman container with a pre-bind receipt, proved actual agent startup removed
  it, verified the running service and exclusive journal, force-killed only that
  exact service's main process, and proved successful post-stop recovery and
  journal reacquisition. No container process was started or image pulled.
- The generated unit passed `systemd-analyze --user verify`, exit 0, without
  installation. Its private artifact remains under
  `/tmp/lcm-managed-agent-BDTsPi/`. The operator guide now documents the commands,
  service installation boundary, exact recovery target and remaining limitations.
- An earlier test run passed its lifecycle assertions but failed while saving
  diagnostics (`chmod` was given an IO handle rather than a filename). That
  harness bug was fixed. The stopped audit-service record was removed only after
  rechecking its invocation, fixed commands, zero main PID and recovered journal.
  Private diagnostics remain in `/tmp/lcm-managed-agent-oEHzgZ/`; subsequent
  successful gates retained `/tmp/lcm-managed-agent-cfWo3H/` and
  `/tmp/lcm-managed-agent-BDTsPi/`. None retains a live audit service or container.
- `p4-managed-host-final.log` exits **2** with `cpu_controller_missing`, as
  required. Final read-only inventory confirms zero audit services and the
  unchanged original 16 containers. `git diff --check` passes.

**P3/P4 remain in progress.** This is real managed-agent death and stopped-resource
recovery evidence, not a running numerical/terminal-container death or effective
quota certification. CPU delegation and an actual Docker Engine remain absent.
Native-unit ownership/recovery, remote preparation/job service and protected
API/client integration, and allowlisted terminal/PTY ownership are next required
work. `assigned_execution` and `private_terminal` remain false. P5 portable
terminal, P6 scientific consumers and P7 release verification retain their full
approved scope; the active goal is not complete.

## P4 — native service identity and recovery checkpoint

`NativeRecovery.jl` now implements the second physical cleanup backend. It uses
the existing bounded command owner and typed systemd inspection, with a fixed
local Unix bus address. `native_scope` verifies socket ownership, bus machine
identity and user-manager ownership before hashing the local machine/user
namespace. Native receipts still use the same private journal, full assignment
fence and immutable binding format as container receipts.

`native_resource_description` derives fixed ownership metadata from the shared
receipt labels. `inspect_native_unit` checks the generated UUID name, transient
marker, invocation, user `app.slice`, main-process state and finite whole-group
shutdown policy. `remove_owned_native!` binds a started crash-gap intent before
control, rechecks identity/scope before mutations, and never signals a saved PID
or clears global jobs/failed services. `recover_native!` attempts every native
record; the actual `recover-agent` dispatcher now reconciles native and container
receipts before startup can announce availability.

Real systemd testing exposed and resolved two important lifecycle cases:

- A queued, never-started service has no invocation yet; the installed manager
  encodes it as an empty byte array. Only an unbound receipt with matching fixed
  metadata, inactive state, zero PID and no assigned control group can use this
  pre-activation path. It cancels only that unit's queued start. A started or
  bound resource cannot fall back to the exception.
- A stopped unit definition may remain cached because another unit references
  it. Recovery does not remove the referencing unit to force garbage collection.
  It requires matched inactive metadata, no process, no pending job and an empty
  or removed exact kernel group on an actual cgroup-v2 filesystem. Repeated
  cleanup is read-only after the receipt is retired. Failed inspection, changed
  identity, pending jobs and surviving descendants keep ownership unresolved.

The Julia-docstring skill was used to keep cleanup identity separate from lease
authority, effective resource limits and model preparation in the new APIs. The
operator guide documents the native cleanup entry points and these boundaries;
native scientific launch remains explicitly unavailable in the current driver.

Verified evidence under `/tmp/lcm-runtime-platform.9uEfiTvj/`:

- `p4-native-recovery-runtime-verified.log`: **1,281** runtime assertions pass,
  exit 0. This includes **79** native scope, ownership, replacement, queued-start,
  cached-record, failure-retention, idempotence and mixed-backend assertions.
- `p4-native-recovery-systemd-verified.log`: **24** real lifecycle assertions
  pass, exit 0. The test used only receipt-owned, finite `sleep` services: it
  removed one while proving another unchanged, exercised the actual agent root
  recovery after journal closure, canceled a queued unit before process start,
  then removed the start-barrier service and its helper group. No scientific or
  terminal executor was admitted, and no image/container was created or started.
  Final private diagnostics are in `/tmp/lcm-native-recovery-Mwm6TG/`.
- Earlier direct/restarted/root-dispatch variants passed **17** checks each;
  their diagnostics remain in `/tmp/lcm-native-recovery-EZJSoH/`,
  `/tmp/lcm-native-recovery-CLhjWQ/`, and `/tmp/lcm-native-recovery-ALAtDa/`.
- The queued-start test initially used a transient job-timeout property that
  this systemd rejects before launch. Only that unsupported fixture property was
  removed; its finite start barrier and 120-second maximum service runtime
  remained. No executor admission limit was relaxed. Its private diagnostic is
  `/tmp/lcm-native-recovery-dMECiU/launch-error.txt`; that trial and
  `/tmp/lcm-native-recovery-fIactu/` left no service behind.
- The first actual pending-invocation trial correctly retained ownership while
  its encoding was unsupported. `p4-native-recovery-retained-cleanup.log` proves
  the exact `/tmp/lcm-native-recovery-U4bgN3/journal` was reconciled after the
  compatibility fix. Its private evidence remains; no live service remains.
- Final service inventory contains no `lcm-exec-*.service` or audit-agent unit.
  The original 16 containers are unchanged. `git diff --check` passes.

**P3/P4 remain in progress.** Native cleanup is implemented and verified, not
quota-enforced scientific launch. Next work must add the native launch policy,
effective child-entry verification and scientific-driver integration, then the
remote preparation/job service and protected API/client path. Shared terminal
ownership/PTY support, P5's portable terminal, P6's scientific consumers and P7's
full release gates remain required. The host still lacks delegated CPU and a
real Docker Engine; those verification gaps remain explicit. The public
`assigned_execution` and `private_terminal` capabilities remain false, and the
full approved goal remains active.

## P4 — shared native/container scientific driver and guarded native entry

`ManagedScientificDriver` now owns both approved scientific backends through the
same journal, capacity accounting, exact-fence handles and required six-hook
interface. Native-only configurations do not probe a container engine.
`ContainerScientificDriver` remains a compatibility alias, not another owner.
Backend dispatch is internal; it does not replace the root lease, preparation,
authorization or cleanup orchestration. The CLI now uses the shared driver.

`NativePolicy` renders a fixed local transient service, bound to the managed
agent's lifetime, with finite CPU/memory/task limits, zero extra swap, private
quota-limited tmpfs scratch, read-only ordinary filesystems, private user/device/
network/IPC namespaces and finite descriptor/core/message-queue limits. Acquisition
records its intent before constructing or starting the attached service command.
The command uses `env -i` and a fixed Julia entry path; it does not inherit broker
credentials or accept browser source, command, mount or environment arguments.
Source/package identity is rechecked against the approved fingerprint.

The fixed pre-import native guard checks actual kernel identity, exact generated
cgroup, quotas, writable mount capacity, privilege flags and network evidence.
`ExecutorLimits` and cgroup quota verification are shared with container entry;
`ContainerLimits` remains an alias. A failed guard exits 78 before scientific
package loading. Successful results additionally require the exact service
invocation, command, manager dependency and policy to pass post-execution checks.
Native remains trusted scientific execution, never arbitrary terminal isolation.
The Julia-docstring conventions keep these ownership, effective-limit and
preparation distinctions explicit. No widget, theme or application contract was
forked for this backend.

Evidence under `/tmp/lcm-runtime-platform.9uEfiTvj/`:

- `p4-native-driver-runtime-stable.log`: **1,357 runtime assertions**, exit 0.
  Includes **62 native-policy** checks (fixed command/credential filtering,
  host prerequisites, invocation binding and policy drift), **79 native recovery**
  checks, native admission failure and real journal acquisition with no process
  start. The acquisition fixture fingerprints the actual installed scientific
  source without importing the scientific package.
- `p4-native-driver-core-stable.log`: **205 execution-core assertions**, exit 0.
  Includes **38 native kernel-policy checks** and **12 real native entry-denial
  checks**, plus existing container guard, framing, cancellation and cache tests.
  No unconfined scientific profile or terminal was entered by these denial tests.
- `p4-native-driver-systemd.log`: **31 real control-agent lifecycle assertions**,
  plus **53 existing managed-unit assertions**, exit 0. Both uniquely owned agents
  (container unavailable and native unavailable) remained operational, were
  forcibly terminated, completed exact post-stop recovery and released their
  journals. No numerical/terminal executor or container was started. Private
  diagnostics remain in `/tmp/lcm-managed-agent-kfp1cp/` and
  `/tmp/lcm-managed-agent-lKNmnB/`.
- `p4-native-policy-unit.log`: **3 inert native service-validation assertions**,
  plus **53 managed-unit assertions**, exit 0. The actual generated directives
  pass `systemd-analyze --user verify` without installation or start. Private
  generated files remain in `/tmp/lcm-native-unit-i0SNd6/`.
- `p4-native-host-check.log`: the actual `lcm runtime check-host --runtime native`
  exits **2**, reporting `cpu_controller_missing`, with no resource creation.
- The first full run, `p4-native-driver-runtime.log`, overlapped an edit to
  fingerprinted core source and failed during a preparation check. The complete
  unchanged-source rerun above passes, including that same lease-loss test. No
  production timeout, source check or isolation requirement was weakened.
- Final read-only inventory contains no executor or audit-agent test service;
  all original 16 containers remain unchanged. `git diff --check` passes.

**P3/P4 remain in progress.** This checkpoint implements guarded native launch;
it does not certify a successful numerical native/container launch on this host,
which still lacks delegated CPU and an actual Docker Engine. Next work is the
remote preparation/job service and protected API/client integration, followed by
shared terminal ownership/PTY support, P5's portable terminal, P6's scientific
consumers and P7's release gates. `assigned_execution` and `private_terminal`
remain false until their actual end-to-end paths exist. The approved goal remains
active, with no required phase removed.

## P3/P4 — remote preparation, expiring evidence and shared controls

`ScientificCoordinator` and `AgentScientificService` now carry explicit
preparation, retained-state inspection and exact-request cancellation through
separate bounded NATS connections. Heartbeat/lease schedulers and durable job
transport remain independent. Both owners are inert at construction and are
started/stopped by the existing root services; they reuse the existing SQLite
assignment authority and `ScientificResources`, not another lease registry.

New strict passive protocol records fence commands/reports to the run, owner,
worker boot, coordinator, profile and assignment generation. Increasing revisions
reject reordering; bounded mutation-ID/digest tombstones prevent an uncertain HTTP
retry from preparing twice. Status queries coalesce, do not create a child or
implicitly prepare, and do not consume mutation history. A fresh ready report
requires actual child inspection. The coordinator subtracts broker latency and
limits readiness to the existing report/lease lifetime; the browser subtracts
HTTP latency and expires its own display without depending on another poll.
An unexpired previous report remains readable during a successor status query.
Mutations clear it immediately. At sub-millisecond remaining validity, the public
state is unknown rather than an inconsistent ready/zero-validity response.

The protected `/runtime/api/assignments/UUID/science` endpoint implements GET
status and explicit POST prepare/cancel under the existing owner, origin and CSRF
checks. Foreign owners are rejected before request-body parsing. Inputs are
passive and limited to 64 KiB; image/command/code execution fields are not an API.
Profile validation and mandatory launch isolation still apply in the agent.
`preparation_control` is advertised separately; `assigned_execution` and
`private_terminal` remain false.

`PreparationStatus(client, role; parameters=Dict())` now offers Prepare executor
and Cancel preparation using the same browser client/renderer as the control
panel and presentation/workbench Bonito views. Polling stays independent of
inventory and diagnostics, is bounded, and is removed on teardown. Hidden pages,
expired/released assignments, stale inventory and pre-mutation responses cannot
retain browser readiness. Background queries do not permanently disable the
Prepare button or appear as cancellable preparation. The existing palette, form
and button CSS remain authoritative; no new theme or page-specific style was
introduced. X-ray exposes action names but excludes preparation input values.
The Julia-docstring skill kept the acceptance, preparation, inspection and
redaction contracts explicit in the component/interface documentation.

Verified evidence under `/tmp/lcm-runtime-platform.9uEfiTvj/`:

- `p4-science-protocol.log`: **133 protocol assertions**, exit 0; **34** new
  command/report bounds, strict shape, passive-input and ready-evidence checks.
- `p4-science-runtime-stable.log`: **1,443 runtime assertions**, exit 0;
  coordinator correlation/retry/expiry and **32** agent preparation/cancellation
  checks use the existing finite test-only child driver. The later millisecond
  boundary refinement passes **242** focused inventory/lease/coordinator checks
  in `p4-science-expiry.log`, including **51** preparation-coordinator checks.
- `p4-science-tls-final.log`: **162 actual TLS assertions**, exit 0, including
  **8** new worker-identity/subject-permission checks and **35** checks that
  prepare, inspect, retry and cancel a finite child through the protected API
  and production root schedulers. Closing only the scientific connection leaves
  heartbeat, lease control and HTTP health responsive. Resource-confirmed release
  removes the child. The first run, `p4-science-tls.log`, also passes 162.
  The finite test-only child is not an effective OS-quota certification.
- `p4-science-http.log`: **94 protected HTTP/offline assertions**, exit 0;
  owner/CSRF/extra-field enforcement and an unavailable scientific channel do not
  allocate a UI host or executor.
- `p4-science-publisher.log`: **1,060 publisher assertions**, exit 0, including
  **36** shared-runtime-component, **447** visual-contract, **40** workbench and
  **19** X-ray preview checks. The earlier focused Bonito wrapper test also passes.
- `runtime_client.mjs` passes with finite expiry, latency subtraction, independent
  polling, lost/pre-mutation responses, hidden-page invalidation and teardown.
  `p4-science-browser-final.log` passes in actual separate headless Chrome with
  light/dark/light controls, 1280/390-width layout, explicit preparation/progress/
  cancellation, ready-to-released transition, uncertain-action retry and detach.
  Screenshots are retained in `/tmp/lcm-runtime-controls.ihuFpF/`; the earlier
  browser run also passes in `/tmp/lcm-runtime-controls.aumJ8S/`.
- Both uniquely named broker fixtures were removed, along with their private test
  credentials/certificates. Non-secret diagnostics remain in
  `/tmp/lcm-runtime-control.K3gAcL1j/` and `/tmp/lcm-runtime-control.pfxpg9oK/`.
  Final inventory contains the same original 16 stopped containers. No managed
  agent service, numerical container or terminal was launched. `git diff --check`
  passes.

**P3/P4 remain in progress.** Next implementation is the production durable job
service: connect the existing `BrokerJobs`/`AssignedDelivery` transport to
`execute_assigned!`, preserve result-before-ack and bounded redelivery, and add
protected submission/result/cancellation plus job/stage correlation to the shared
client. Do not redo the now-implemented remote preparation channel, native driver
or resource recovery. Shared terminal ownership/PTY support, P5's portable Julia
terminal, P6's real scientific deck/workbench and P7's release gates remain
required. This host still lacks delegated CPU and an actual Docker Engine; those
external verification gaps remain explicit, not reasons to weaken isolation or
remove an approved phase. The full goal remains active.

## P3/P4 — exact prepared targets, durable agent execution and job receipts

Assigned v2 jobs/results now carry a strict `PreparedExecution` record containing
the actual executor UUID, process generation and preparation key. The agent
checks that target under the same lock as admission; delayed work cannot run
against a replacement child or a different prepared model. The earlier,
unreleased v2 message shape is rejected rather than guessed. Legacy v1 remains
unchanged. The child returns its registered result-schema version through the
shared framed executor. `ScientificOutput` retains that schema, normalized value
and warnings together, independently of later status inspections. Provenance
explicitly labels the approved environment digest rather than inventing an
engine package version.

`AgentJobService` now connects the existing worker-specific durable consumer to
`ScientificResources`. Its connection, bounded polling and per-delivery tasks
are independent of heartbeat/lease and preparation traffic. Construction remains
inert; the root agent starts it only after resource recovery. Delivery payload,
subject, acknowledgement metadata, worker identity and prepared target are
checked before numerical admission. Capacity includes results awaiting storage.

The service reads the exact durable result before executing. A matching result
is acknowledged without another child request. A redelivery with no stored
result returns `execution_uncertain`, never an automatic numerical rerun.
Progress acknowledgements extend the delivery timer while authority is live;
terminal acknowledgement follows result persistence. Persistence failures retry
the retained outcome, not the calculation, within the live lease and five seconds
beyond the job deadline. Independent shutdown cancels only its exact pending
requests and joins owned tasks; the shared resource owner still performs physical
lease cleanup. Inline results are limited to 256 KiB; larger values currently
produce `result_payload_limit`. Assigned artifact delivery remains required work.

Real TLS testing exposed a race between background readiness inspection and
job admission. The shared owner now waits at most five seconds behind a
read-only inspection and atomically rechecks the target before admission. It
does not queue behind another operation or preparation. A deterministic test
holds a real inspection in progress while a delivered job waits for admission.

SQLite schema **4** adds `JobRecord` submission receipts, not a second result
store or lease authority. Transactions preserve the original job ID, absolute
deadline, prepared target and inputs for an identical request ID, including
after restart/release. Changed-input reuse and cross-owner reads fail. Admission
allows one pending job per lease and 256 receipts per assignment lifetime;
history reads are bounded. Explicit migration from schemas 1/2/3 retains the
existing run/worker/lease bookkeeping and creates a private consistent backup.
No operator database was migrated. `runtime/CONFIGURATION.md` documents these
contracts and remaining limitations. The Julia-docstring skill kept the new
records, return values, lifecycle and low-level authority caveats explicit.

Stable-source verification under `/tmp/lcm-durable-jobs.0oTwUtZj/`, all exit 0:

- `protocol.log`: **140 assertions**, including strict prepared-target shape,
  identity/generation/hash validation and assigned result provenance.
- `core.log`: **205 assertions**; `worker.log`: **80 assertions**, including
  existing supervised cancellation, framed output and real engine parity.
- `runtime-stable.log`: **1,529 assertions**, including **37** inert job-owner/
  result-contract checks, **42** transactional receipt/history/migration checks,
  and the scientific owner tests with exact target rejection and child-reported
  schema 1.2. Earlier `runtime-receipts.log` also passes 1,529 before the final
  independent job-service shutdown refinement.
- `tls-stable.log`: **241 actual TLS assertions**. These include **52** protected
  API/root-scheduler preparation, execution, cancellation and control-independence
  checks, and **62** durable-agent checks for malformed deliveries, lost input
  acknowledgement recovery, uncertain execution, stale targets, inspection
  contention, a 32-second operation across the provisioned 30-second AckWait,
  and independent shutdown without a surviving child. Broker diagnostics remain
  in `/tmp/lcm-runtime-control.xbIuLN2X/`.
- `publisher.log`: **1,060 assertions**, including **447** visual contracts,
  **40** workbench and **19** X-ray preview checks. No CSS or widget-specific
  rendering variant was introduced in this checkpoint.
- Earlier TLS failures are retained in the same log directory: the real
  inspection/admission race was fixed; the manually controlled long-delivery
  fixture was corrected to use real monotonic time and independent coordinator
  presence updates. The final run used stable source and passes shutdown too.
- `git diff --check` passes. Every uniquely named test broker and its private
  credentials/certificates was removed. Final Podman inventory is the same
  original **16 stopped containers**. No managed numerical/terminal container,
  installed agent service or public deployment was started. Finite test-only
  children do not certify this host's missing CPU-controller enforcement.

**P3/P4 remain in progress; the full goal remains active.** The durable agent
path and receipt store are implemented, but the coordinator has not yet connected
those receipts to protected submission/result endpoints and the shared browser
job client. `assigned_execution` and `private_terminal` therefore remain false.
Next: add the bounded coordinator job publisher/result reconciler using these
receipts, protect result retrieval by owned receipt plus exact provenance, and
complete cancellation for both queued and running jobs. A cancellation that
arrives before delivery must prevent a later start; merely forwarding the current
request cancellation is insufficient for that queued case. Add job/executor/stage
diagnostic correlation and the existing artifact path, then shared client/job
controls and last-good-data provenance. Do not redo remote preparation, native
drivers, the durable worker consumer or schema-4 receipts. Shared terminal
ownership/PTY support, P5's portable terminal, P6's scientific consumers and P7's
release gates remain required. Actual Docker and effective numerical/terminal
CPU-limit verification remain explicit external gaps, not reasons to weaken the
approved isolation contract or remove phases.

## P3/P4 — protected job orchestration and durable queued cancellation

The production `ControlService` now owns a separate `JobCoordinator` using the
existing SQLite receipts, lease coordinator, preparation channel and diagnostic
buffer. Submission saves the exact approved request before asynchronous delivery;
the HTTP caller never waits for calculation. Four fair, finite reconciliation
tasks and four independently bounded HTTP result reads cannot occupy heartbeat
or preparation scheduling. Shutdown drains both kinds of task before closing the
job connection. Identical explicit submission retries return the original receipt
without extending its deadline or choosing a different executor.

The gateway now provides protected submission, run job history, receipt, result
and cancellation endpoints. Ownership is checked before input parsing or broker
access. Durable result acceptance verifies the complete saved fence, prepared
target, job UUID, operation and input hash. Broker absence is distinct from an
actual missing result. Lost authority advances bookkeeping even while offline;
an unresolved deadline becomes uncertain, not successful. Later durable evidence
may resolve uncertainty but can never requeue the operation.

Cancellation now covers both running and not-yet-delivered jobs. SQLite schema 5
adds owned cancellation intent and a distinct worker-acknowledged flag. Workers
retain bounded exact-assignment cancellation tombstones across child replacement
and re-preparation. Cancellation-before-publication first records that tombstone,
then delivers the original request for an actual durable canceled outcome. It
does not fabricate a result from a control ACK. Saved completed results prevail
over late cancellation. Only this idempotent `cancel_job` action may be safely
resent after a lost/rejected reply; uncertain preparation remains non-replaying.

Migration from schemas 1–4 stays explicit, stopped-coordinator-only, backed up and
transactional. No operator database was migrated. New job/executor/stage fields
use the shared bounded owner-filtered diagnostics and renderer. Scientific input,
output, credentials and arbitrary exception text are excluded.

Shared browser `RuntimeClient` now supplies `submitJob`, `listJobs`, `job`,
`jobResult` and `cancelJob`. It rejects non-passive/unbounded inputs, mismatched
receipt/result identities and changed result executor generations. Invalid
mutation replies retain their original retry identity and remain uncertain.
Read APIs never resubmit or follow worker-supplied artifact URLs. These methods
are reusable plumbing, not yet the scientific execution-control/view component.
No additional theme tokens or page-specific styles were introduced.

Verification under `/tmp/lcm-job-api.YCAXUIJ6/`, all exit 0:

- `import.log`: inert runtime import, SQLite schema 5.
- `protocol.log`: **143 assertions** including queued-cancellation wire shape.
- `runtime.log`: **1,627 assertions**, including **30** worker cancellation,
  **19** cancellation persistence/migration, **39** coordinator ownership/drain/
  uncertainty checks and **10** cancellation retry-ordering assertions.
- `tls.log`: **274 assertions**, including **85** protected API/production
  scheduler checks. Tests submit through real HTTP, reconcile through authenticated
  NATS, execute finite fixture children, retrieve exact results, cancel a running
  operation and deterministically cancel before publication. Cross-owner reads
  and malformed writes are denied; result transport loss leaves health/control
  usable. Prior durable redelivery/lost-ACK and long-operation progress tests pass.
- `client.log`: shared browser client regression, including lost mutation reply,
  explicit immutable retry, invalid receipt identity, wrong executor generation,
  history after revocation and finite passive inputs.
- `controls-browser.log`: real isolated Chrome shared-renderer gate passes draft
  retention, offline pinning, both themes, compact layout, explicit retry,
  registration, release, diagnostics and teardown. Browser diagnostics:
  `/tmp/lcm-runtime-controls.3yHtDV/`.
- `publisher.log`: **1,060 assertions**, including **447** visual contracts,
  **40** workbench checks and **19** X-ray preview checks.
- `git diff --check` passes. The temporary NATS-only container and its private
  certificates/password files were removed. Remaining transport diagnostics are
  `/tmp/lcm-runtime-control.xK5Z58yr/` (`nats.conf` and `broker.log` only). Podman
  still lists precisely the original **16 stopped containers**. Test cancellation
  intentionally terminates finite Julia fixture children; this is not a quota-
  enforced numerical/terminal launch or an actual Docker certification.

**P3/P4 and the full P0–P7 goal remain active.** Next is assigned artifact delivery
and the shared scientific execution/view controls with last-good provenance.
Reuse the existing content-addressed artifact infrastructure; private runtime
results must be accessible only through an owned job receipt, never the legacy
public `/artifacts/sha256` route. Remote agents must not require shared local
storage. Current oversized inline assigned results still return an explicit
failure. `assigned_execution` and `private_terminal` remain false. Then complete
terminal container/PTY ownership, P5's shared terminal, P6's real deck/workbench
consumers and P7's aggregate release verification. Do not redo the now-verified
job API, cancel tombstones or migrate an operator database implicitly.

## P3/P4 — private assigned artifacts and remote storage

The production agent now stores successful scientific values above 64 KiB using
an optional private artifact backend. The coordinator retrieves them only through
an owned job receipt and matching durable result. `GET`/`HEAD
`/runtime/api/jobs/UUID/artifact` verifies metadata, length and SHA-256 before
serving JSON with no-store, nosniff and attachment disposition. Ownership precedes
broker/storage access. The legacy public digest route remains unchanged and does
not expose this private namespace.

Shared S3 configuration, key construction and metadata definitions were extracted
from the existing publisher/worker copies into `common/artifact_contract.jl`.
Those APIs and their existing public layouts remain compatible. Runtime source
declarations and the legacy worker fingerprint include the extracted file.
Existing deployment copies include it; scientific profile images and the
execution core still neither import the storage SDK nor receive credentials.
Runtime uses the existing pinned AWS SDK for signing, with the installed HTTP
client providing bounded transport instead of implicit SDK retries.

Both server-owned control and agent configurations accept the same optional
`artifacts` table. Local mode uses a dedicated private filesystem root. Remote
mode uses verified HTTPS S3 storage, private credential files and job-specific
worker/run/boot/lease/generation prefixes, without a shared artifact filesystem.
Agents retain write/delete-only permissions for their own worker prefix; the
coordinator has read permission, not write/delete. Anonymous access is denied.
Neither browser input nor a worker-supplied retrieval URL selects a storage path.

Storage requests disable proxy inheritance, redirects, cookies and automatic
retries. Each has a five-second deadline and bounded response buffer. Four
independent coordinator artifact reads are joined during shutdown without holding
heartbeat/status scheduling. `RuntimeClient.jobArtifact` uses the private job URL
and its own finite 15-second read, leaving inventory polling unchanged. No new
palette, theme variant, page-specific renderer or X-ray value exposure was added.

Uploads commit metadata last. Normal failures attempt cleanup of only the job's
two exact objects and produce a durable `artifact_unavailable` failure without
rerunning computation or falsely invalidating a still-prepared child. Filesystem
temporary writes are private, synchronized and atomically renamed. Remote uploads
have no local disk staging. SDK, TOML and transport exceptions are discarded before
raising the public storage error, including nested failed-task diagnostics.

Bounds and limitations are explicit in `runtime/CONFIGURATION.md`: storage and
downloads accept at most 4 MiB; the current executor still has its tighter 1 MiB
complete framed-response limit. Successful artifacts intentionally survive lease
release. An unavailable remote deletion may leave an orphan until operator-owned
retention expires; the fixture uses two days, not a production retention default.
Filesystem crash residuals require stopped-owner maintenance. This is bounded
cleanup with an honest failure policy, not a promise of atomic transactions across
SQLite, JetStream and S3.

Verification under `/tmp/lcm-private-artifacts.AOGl316w/`, final runs exit 0:

- `runtime-final.log`: **1,690 assertions**, including **51** private artifact
  configuration, filesystem isolation, corruption, capped transport, failed-write
  cleanup and nested credential/endpoint-redaction checks. Reader shutdown and
  both real configuration parsers are also covered.
- `tls-final.log`: **306 assertions** against disposable authenticated TLS NATS
  and MinIO containers. Includes **9** actual storage IAM checks and **108**
  protected API/production-scheduler checks: a 60,000-value result exceeds the
  broker's message limit, travels through private storage and is downloaded by
  its owner. Foreign/anonymous access, another worker's prefix, a read-only
  writer, private data through a public digest URL and untrusted TLS are rejected.
  Misconfigured storage yields a durable failure while preserving executor
  identity; subsequent correctly configured work succeeds.
- `protocol-1.log`: **143 assertions**; `core-1.log`: **205 assertions**, including
  its no-AWS/no-engine import boundary and supervised lifecycle tests.
- `worker-2.log`: **80 assertions**, including legacy artifacts and real
  LineCableModels scientific parity. `publisher-1.log`: **1,060 assertions**,
  including **447** visual, **40** workbench and **19** X-ray preview checks.
- `client-2.log`: shared browser client passes private same-origin artifact reads,
  input/path validation, result provenance, finite mutation handling and teardown.
  `browser-1.log`: isolated Chrome shared runtime controls pass in both themes;
  diagnostics `/tmp/lcm-runtime-controls.au4zad/`.
- Earlier focused failures remain as evidence: corrected explicit HTTP client
  proxy configuration, HTTP fixture body conversion, and storage-fixture startup
  ordering. Final stable-source regressions include the redaction refinement.
- `git diff --check` passes. All newly created broker/storage/initialization
  containers, passwords and TLS keys were removed. Final diagnostics remain in
  `/tmp/lcm-runtime-control.oj4tmgiL/`; only non-secret policies/configuration and
  service logs remain. The original **16 stopped containers** are unchanged.
  No operator database migration, installed service, public deployment or managed
  scientific/terminal container was started. Finite test children do not certify
  effective CPU limits; actual Docker verification remains unavailable here.

**P3/P4 and the complete P0–P7 goal remain active.** Private artifact delivery is
now implemented and verified. Next: compose shared scientific execution/view
controls on the existing client, with explicit submission/cancellation, unchanged
retry identity, current-draft/assignment/executor fencing and visible last-good
data/provenance. Then complete terminal process/PTY ownership and relay, P5's
portable terminal, P6's scientific deck/workbench and P7 release verification.
`assigned_execution` and `private_terminal` remain false until their complete
consumer paths are ready. Do not redo job transport, cancellation, artifacts or
the source/physical-supervision foundations, and do not weaken mandatory limits.

## P3/P4 — shared explicit calculations and retained scientific views

Added the exported `ScientificJob` component and `ScientificResult` binding in
`src/widgets/RuntimeControls.jl`, with one `RuntimeJob` intent tracker on the
existing `assets/runtime-client.js` and the existing shared renderer. The
inventory-only gallery demonstrates the same component without allocating a
worker or submitting work. No CSS palette, broker client, result transport or
per-page execution implementation was added.

The control captures explicit Run intent, displays durable state, separates
cancellation request/acknowledgement from completion, and retries unconfirmed
mutations only with the original request UUID and inputs. Input changes never
submit work. It retains one last-successful value/provenance, marks it outdated
when inputs or evidence cease to match, and rejects late completion after a
known input/assignment/executor replacement or revocation. Expiring readiness
cannot create a null target at click time. Unknown evidence and an offline broker
cannot establish currentness or enable Run.

The Bonito bridge validates passive result shape, strict assigned provenance,
run/role/operation, current input hash and a render-local draft token. Invalid
inputs are represented explicitly by a nullable wire value, not an Observable
conversion exception. Identical input echoes retain the same draft token.
Private inputs/results stay outside X-ray and ordinary `show` output. This is
display state, not an authorization mechanism. Gateway enforcement is unchanged.

Shared-client ownership now also covers departing frames: non-BFCache `pagehide`
closes their controls, while removing one ordinary component preserves sibling
subscriptions. Invalid execution configuration neither alters the DOM nor leaves
an unused client registered. No teardown implicitly cancels accepted computation.

Verification on 2026-09-08, diagnostics `/tmp/lcm-execution-controls.OfVFQfl4/`:

- `publisher-1.log`: **1,101 assertions passed**, including **77** runtime
  component/result-projection checks, **447** visual, **40** workbench and **19**
  X-ray checks. `julia-controls-2.log` independently passes the focused 77 checks.
- `jobs-final.log`: Node tests exercise the real shared client/tracker against an
  explicit HTTP fixture: inert construction, input changes during private artifact
  download, exact retry identity, cancellation, executor replacement, revocation,
  readiness expiry, stale evidence and sibling-safe cleanup.
- `client-final.log`: existing shared-client regression passes. The final tracker
  change does not alter its transport implementation.
- `controls-browser-1.log`: existing isolated Chrome control-page test passes in
  both themes; diagnostics `/tmp/lcm-runtime-controls.DmgljU/`.
- `bonito-3.log`: actual gateway/Bonito sessions and binary sockets pass the
  scientific input/result round trip, persistent output identity, last-good
  retention, stale-result fencing, explicit retry, invalid drafts, light/dark
  gallery parity and removed-frame cleanup. Its scientific HTTP transport is
  **explicitly mocked** in a test-only client; this is not a new broker/solver
  execution claim. The previously verified production API/worker/artifact tests
  remain the separate authority for that transport boundary. Existing same-run
  frames, independent UI hosts, reconnect, X-ray, stop/survival and clean restart
  checks also pass. Diagnostics `/tmp/lcm-runtime-bonito.wqGlznqB/`.
- Earlier fixture failures remain recorded: incorrect lease/run field order,
  top-level-await syntax and unavailable global theme helper. The nullable input
  binding failure was an implementation defect, now covered by Julia and browser
  regressions. All isolated gateway, UI-host and browser processes were joined by
  their test owners. No container, operator database, service or deployment was
  changed. `git diff --check` passes.

**The complete P0–P7 goal remains active; this checkpoint is verified progress.**
The shared calculation/view control is now available. Next required work is
allowlisted terminal-container/PTY ownership and relay (P4/P5), then the concrete
scientific deck/workbench (P6) and aggregate release verification (P7). Keep public
`assigned_execution` and `private_terminal` false until their complete consumer
paths are ready. Reuse the existing physical resource owner and private transport;
do not rebuild job/artifact/control foundations or weaken unavailable CPU limits.

## P4 — bounded local PTY transport and real Julia REPL behavior

Implemented internal `runtime/src/TerminalProcess.jl`. This is the actual PTY
transport needed by the terminal driver, **not** a native-terminal shortcut or a
completed private terminal. `TerminalProcess` construction is passive; a retained
handle precedes `start_terminal!` and every partial acquisition. Production
callers still need the approved container driver, current lease and physical
resource journal. No browser/CLI terminal entry point or public capability was
enabled, and no second command parser or Julia evaluator was introduced.

The implementation uses a nonblocking Linux PTY master, bounded FileWatching
polls and one I/O pump per terminal. It does not use libuv's potentially blocking
PTY-master write fallback. It retains output in a fixed-size byte ring, reports
cursor gaps explicitly, bounds input queue/chunk sizes, and enforces output-rate,
stalled-write and command-lifetime deadlines. Raw byte boundaries preserve UTF-8
and escape-sequence streaming. Metadata display omits command/input/output data.
Resizing updates terminal cells and signals the original attached process.

Natural exit, input/output failure, lifetime expiry and explicit close join the
original command and I/O watcher before closing the descriptor. Incomplete
cleanup retains the same handle/descriptor for retry. Input chunks are cleared
when written or retired. Attached CLI retirement does not claim container deletion
or lease release. The only new runtime dependency is Julia's already-installed
FileWatching standard library; its environment was resolved offline.

Verification on 2026-09-08, diagnostics `/tmp/lcm-terminal-pty.RsIUiuJ2/`:

- `pty-3.log`: **112 assertions passed** with two Julia threads. Actual PTYs and
  finite Julia children verify Unicode, multiline definitions, completion,
  session history, window-size changes, Ctrl-C recovery, distinct REPL state,
  explicit exit, output floods, stalled input, lifetime limits, partial cleanup
  retry and repeated descriptor retirement without growth.
- `pty-single-thread.log`: the same **112 assertions pass with one Julia thread**,
  including a surviving sibling while another terminal floods or cannot read.
- `runtime-1.log`: the full runtime suite passes **1,802 assertions**. Its existing
  scientific/control/identity/recovery tests remain intact; subsequent source
  edits were docstrings only. The PTY suite is now included in `test/runtests.jl`.
- The REPL fixture is explicitly test-only and has an independent 60-second
  kernel alarm. It uses Julia's actual REPL and interactive SIGINT behavior.
  The first run's script-mode SIGINT exit was corrected in the fixture; no
  transport workaround masks it. No numerical/terminal container was started.
- `host-check.log`: the existing CLI still fails closed with
  `cpu_controller_missing` for rootless Podman. `containers.log` confirms all
  **16 original containers remain stopped and unchanged**. Actual Docker and
  running container limit/resize propagation are not verified on this host.
- All finite fixture process/PTY owners completed cleanup. No service, operator
  database, broker policy or public deployment was changed. `git diff --check`
  passes.

**P4 and the complete P0–P7 goal remain active; this turn made verified progress.**
Next: integrate this transport with terminal lease/container ownership, reusing
the same managed physical journal and total capacity used by scientific profiles.
Do not open a competing terminal journal or duplicate AgentService. Current
`ManagedScientificDriver` advertises only scientific profiles and explicitly
marks terminal profiles unavailable; `ScientificResources` also requires a
scientific-only registry. Integrating both kinds needs shared physical ownership
and kind-specific resource dispatch, followed by sole-writer session authority,
disconnect/idle expiry, separately bounded ordered relay and the portable xterm
component (P5). Then complete P6 consumers and P7 acceptance. Public
`assigned_execution` and `private_terminal` stay false until complete consumer
paths are verified; missing CPU delegation is not permission for native REPL fallback.

## P4 — shared managed ownership and lease-bound terminal sessions

Scientific and terminal resources now compose through one physical owner and
one agent ledger. `ManagedResourceDriver` is a compatibility alias for the
existing managed driver, not a competing supervisor. Its per-lease handle retains
either the existing scientific supervisor or a terminal PTY. Both kinds use the
same capacity limit, durable journal, fixed container policy and exact receipt
cleanup. `ManagedScientificView` and `ManagedTerminalView` expose only their
verified kind-specific profiles and release only their own handles.

Borrowed views close admission and join in-flight acquisition before collecting
cleanup targets. The admission condition does not serialize unrelated physical
operations; a partition waiting for one acquisition does not block the other
partition. Closing a borrowed view does not close its sibling or the parent
journal/runner. `ManagedAgentResources` composes both logical owners, binds both
to the existing `AgentLeaseLedger`, dispatches release by profile kind, and closes
both partitions and the parent with retryable failure. The existing CLI now uses
this root. `AgentService` obtains its scientific partition by dispatch; there is
no second agent scheduler or terminal lease table.

`TerminalResources` now provides the internal session boundary:

- Passive construction/binding; explicit open under a usable exact assignment.
  One fresh stream UUID and one writer UUID; matching retries do not create a
  second process. A different writer cannot take over a disconnected stream.
- Fixed container attachment by full owned ID, after stopped acquisition and
  physical policy verification. The successful container entry guard installs
  a receipt-bound Julia REPL initialization marker. The session owner observes
  that exact bounded marker, rechecks physical policy/process liveness and strips
  the startup prefix before admitting input. No user-code evaluation occurs in
  the agent and no prompt-string heuristic establishes readiness.
- Consecutive bounded input chunks. An exact retry of the last admitted sequence
  is acknowledged without replay; altered/older/skipped input fails. Lease
  revocation is atomic with local queue admission. An ACK is not a claim that
  Julia evaluated the input. Ctrl-C uses the same sequenced byte path.
- Bounded relative output cursors with explicit gaps; status contains only
  bounded state/counters and sanitized failure codes. Commands, writer tokens,
  input and output are not ordinary diagnostics or X-ray metadata.
- Finite startup, disconnect, idle and cleanup deadlines. Disconnect/idle expiry
  applies during startup as well as after readiness. Passive reads and repeated
  disconnects do not renew activity. Only the same writer can reconnect within
  grace; expired/exited sessions never implicitly restart.
- Lease loss, startup failure, natural exit and timeout stop the original PTY and
  attempt exact physical retirement. A refused container removal does not leave
  terminal transport running. Failed cleanup retains ownership for retry; ended
  session identities remain until lease release.

The shared image recipe includes `terminal-ready.jl`. Rebuilding and approving a
new image digest is still an operator action, not performed by this change.
Public `assigned_execution` and `private_terminal` remain false. Terminal
profiles remain excluded from production advertisement with the more accurate
`terminal_relay_unavailable` reason: the separately authorized relay and browser
component are not yet connected. No native terminal fallback was introduced.

Verification on 2026-09-08, diagnostics `/tmp/lcm-terminal-owner.4gQdilRT/`:

- `partitions-2.log`: focused shared-capacity, cross-kind rejection, cleanup,
  shutdown/acquisition-race and fixed-attachment checks pass. Subsequent full-suite
  coverage adds the composed-agent/ledger/scientific-service assertions.
- `sessions-1.log`: the initial **94 terminal session assertions** pass with two
  Julia threads, using the actual `--interactive` Julia REPL. They demonstrate
  two independent Julia namespaces, excluded second writers, non-replayed input,
  Ctrl-C recovery, reconnect, deadlines, lease expiry and cleanup retry.
- `sessions-single-thread.log`: the finalized **106 terminal assertions** pass
  with one Julia thread, plus the existing 55 lease assertions. This includes
  exact/partial/oversized startup marker checks and disconnect expiry while
  startup is still pending.
- `runtime-1.log`: the first integrated runtime pass completed **1,970 checks**.
  After the final startup-expiry/input-admission tightening, `runtime-2.log`
  completed **1,981 assertions**, exit 0. It includes the existing scientific,
  job, control, identity, recovery, dependency and image-layout regressions, as
  well as the 112 earlier PTY transport assertions. Expected scientific-child
  termination traces from cancellation tests are not terminal-session failures.
- Early failures were corrected and remain in their diagnostic files: an
  unqualified RequiredInterfaces helper in the ad-hoc import command; a sandbox
  precompile-cache permission failure; and an incorrect synthetic container
  constructor/receipt fixture. Final imports, required-interface checks and both
  focused/full regressions pass. No production isolation check was bypassed.
- `host-check.log`: the existing CLI still exits 2 with
  `cpu_controller_missing` for rootless Podman. `containers.log` confirms all
  **16 original containers remain stopped**. No container was created or started.
  Actual Docker and running container PTY/resize/isolation propagation remain
  unverified; finite native test REPLs are not evidence for those guarantees.
- All test-owned Julia processes and PTYs were joined/retired. The two pre-existing
  long-running Julia processes were not touched. No operator database, broker
  policy, existing service or public deployment was changed. `git diff --check`
  passes. Operator documentation describes the implemented and pending boundaries.

**The complete P0–P7 goal remains active; this checkpoint is verified progress.**
Next required work is the separately authorized, bounded terminal relay and
gateway/session integration, followed by the shared locally served xterm
component (P5). Reuse `TerminalResources`, `ManagedAgentResources`, the existing
agent/lease authority and gateway origin/principal checks; do not rebuild PTYs,
scientific dispatch or physical recovery. Enable terminal profile eligibility
only with that complete private path. Then finish the real scientific deck and
workbench consumers (P6) and aggregate acceptance/operator handoff (P7). Missing
host CPU delegation is an external verification limitation, not permission to
weaken isolation or replace terminal containers with native processes.

## P5 — private terminal protocol, broker and gateway relay

Implemented on 2026-09-08 after the shared-owner checkpoint:

- Added transient `TerminalCommand`/`TerminalReport`, exact fenced broker subjects,
  strict 8 KiB byte chunks and 64 KiB frames. Broker grants are worker/direction
  scoped; exact subscriptions and full payload/lease checks preserve assignment
  identity. Private bytes never use JetStream or scientific job records.
- Added separate `BrokerTerminal`, `AgentTerminalService` and `TerminalCoordinator`
  ownership. Root agent/control services compose these with their existing lease
  ledgers. Latest-only command revisions/digests and explicit input sequences
  reject changed retries, duplicates and late replies without reevaluating input.
  Mutation sends are never automatically replayed. Completed/timed-out flights
  release raw command payloads; bounded reports are volatile.
- Added explicit exact-stream stop/restart. Restart joins physical cleanup before
  creating a new UUID/Julia namespace; repeated exact requests cannot create two
  replacements. The existing resource owner remains authoritative.
- Added writer presence separate from idle activity: fresh keepalive every five
  seconds, default 15-second absence detection, then original 30-second reconnect
  grace. Cached keepalive replies do not renew presence. Channel loss/disconnect
  begins grace once; status/output reads do not refresh it. Only same-writer open
  can reconnect before expiry.
- Added the protected assignment terminal WebSocket route. Actual ownership is
  checked before upgrade and every action/output; inventory administrators cannot
  read another owner's private REPL. The browser supplies no fence, revision,
  environment, container command or worker override. One attachment per lease,
  one in-flight action and one waiting action are enforced.
- The HTTP 2.6.6 compatibility adapter bounds its otherwise-unbounded receive
  queue before the reader starts. Four stored frames plus one blocked reader
  frame are bounded; fragmentation/compression are disabled. Socket writes,
  initial attachment, request rate and silence have finite bounds. Abrupt EOF
  closes immediately, not only after a completed WebSocket close handshake.
  Flood cleanup releases the queue and original transport before joining the
  reader. No global HTTP method, installed package or unrelated socket changes.
- Pinned local renderer assets are generated from xterm.js 6.0.0, fit 0.11.0 and
  esbuild 0.25.9. The upstream MIT notices match the installed packages byte for
  byte. `assets/vendor/README.md` records the build; no runtime CDN or npm is
  needed. This is only the renderer dependency, not the finished browser widget.

Verification: `/tmp/lcm-terminal-relay.9ZPDrXWK/`.

- `protocol-2.log`: **198 protocol assertions**, exit 0, including strict
  keepalive grammar, byte/field bounds and redaction.
- `coordinator-1.log`: **225 assignment/lease/terminal coordinator assertions**,
  exit 0; terminal-specific groups contain 58 checks. Owner/admin exclusion,
  exact explicit retries, late-report fences, stream/sequence/cursor consistency,
  private-payload release and shutdown admission are covered.
- `agent-2.log`: **256 lease/session/agent assertions**, exit 0, including 24
  writer-presence checks against finite real REPL processes.
- `runtime-1.log`: the complete runtime test file reached its final test group
  with **2,154 passing assertions** and no failed/error group. The tool session
  was retired during user continuation, so its separate process exit result was
  not retained. It is evidence for those assertions, not the final P7 release run.
- `tls-8.log`: **426 assertions**, exit 0. Includes the direct TLS role/fence
  tests, actual private REPL, bounded-socket flood test and **90 complete
  gateway/control/agent/TLS checks**. The latter verifies two owners' independent
  Julia state, denied foreign/admin/wrong-origin attachment, no duplicate input,
  actual PTY dimensions, stateful reconnect, clean restart, stopped-run closure,
  abrupt TCP loss, unaffected second terminal/public health and no scientific
  job/event contamination. The driver is explicitly a finite native fixture;
  it does not assert container isolation.
- Early fixture failures are retained in `tls-1.log` through `tls-7.log` and
  `sockets-1.log`. Corrections include an empty-set element type, counting a
  blocked Julia Channel put separately from stored entries, fixture keepalives,
  run-scoped job lookup, and EOF-vs-full-WebSocket-close semantics. Cold combined-
  process compilation can consume a finite handshake/lease ACK budget: setup now
  retries connection attempts finitely and confirms exact expired-grant cleanup
  before a fresh reservation. No timeout, CPU requirement or authority check was
  weakened. Cold production startup latency still belongs in P7 acceptance.
- Generated bundle syntax passes `node --check`. SHA-256:
  `058fe53ea7901a1166a2ba6a16a0557e236923457d9e02e7bfbfac21eda072b3` (JS),
  `9a94b33948fa65578113fcc20fa11350c9d0f62272490039c0595beee03730fe` (CSS).
- The TLS harness removed each of its own temporary NATS containers and fixture
  credentials. The most recent diagnostics are
  `/tmp/lcm-runtime-control.2Qs0Uc3F/`. Existing operator containers were untouched.
  Public `private_terminal`/`assigned_execution` remain false; terminal profile
  advertisement remains disabled until the shared component and complete consumer
  gates pass. Missing delegated Podman CPU and actual Docker remain external
  verification limitations, not permission to use native terminal fallback.

**The P0–P7 goal remains active. This is verified progress, not completion.**
Next: implement one shared `JuliaTerminal` browser/Bonito component around the
existing RuntimeClient/WorkerSelector, locally vendored renderer and this exact
gateway path. Exercise both themes, resize, focus/slide-key ownership, output
backpressure, uncertain-input UX and X-ray redaction in gallery/deck/workbench.
Then finish real scientific consumers and full P7 acceptance; do not replace
the already verified resource, lease, broker or private relay owners.

## P5 — shared terminal component and real browser/REPL round trip

Implemented and verified on 2026-09-08:

- Added the passive `JuliaTerminal(client, role; title, rows)` component, one
  browser transport and one locally served xterm renderer. Gallery, live frames
  and workbench composition use those exact assets. The public gallery factory
  defaults to inventory-only; an owned consumer passes its run context explicitly.
  `/widgets/julia-terminal` is documented and present in the rebuilt publication.
- Connect is explicit. Startup, running REPL, disconnection, uncertain action and
  cleanup are distinct. Stop/restart require confirmation; interrupt sends Ctrl-C.
  Clear view changes scrollback, not Julia state. Reconnect preserves only the
  same in-memory writer/stream within worker grace; reload cannot claim recovery.
- One in-flight browser request, 32 KiB unsent input, 8 KiB chunks, bounded socket
  frames and a three-second renderer deadline enforce backpressure. Oversized
  paste is rejected whole. Unmount, departure, stale inventory and changed lease
  close the socket and discard unsent bytes. A lost reply is never replayed;
  reconnect allows inspection, with explicit Resume input or clean restart.
- Readiness text does not mutate unchanged live-region content on every poll.
  Rejected-paste and output-gap notices survive background status updates.
  Delayed connection completion does not steal focus from another control.
- Themes come from the existing brand palette; controls use the shared button
  contract. The fitting utility now receives a border/padding-free inner content
  box. A live screenshot exposed the original last-line clipping despite correct
  cell-count checks; new browser assertions verify the actual rendered boundary.
  Both themes and narrow/wide viewports pass. Terminal keys do not bubble to
  presentation shortcuts; output-controlled title/link/clipboard/palette changes
  are blocked. Selection/copy remains a terminal interaction.
- X-ray exposes public role/title/row configuration, owned CSS and actual action
  names only. No terminal Observable, input/output, writer identity or history is
  introduced into Bonito or its diagnostic metadata.
- The shared runtime client/control panel now distinguishes terminal profiles
  from scientific preparation. It neither polls their scientific endpoint nor
  offers scientific preparation actions. Tests cover mixed role handling and
  retain all existing scientific job/client behavior.
- The terminal eligibility reason is now `terminal_acceptance_pending`, not the
  obsolete `terminal_relay_unavailable`. Public `private_terminal` and
  `assigned_execution` remain false. No container requirement or isolation
  preflight was weakened, and no native terminal fallback was enabled.

Verification and retained evidence:

- `node test/integration/runtime_terminal.mjs`: **47 transport checks**, exit 0,
  using an explicit socket fixture (no broker/REPL claim).
- `node test/integration/runtime_terminal_browser.mjs`: **39 browser checks**,
  exit 0, actual Chrome/xterm plus explicit HTTP/socket fixtures. Latest screenshots
  are `/tmp/lcm-terminal-browser.W0e6Sq/`. Covers both themes, frame bounds,
  keyboard ownership, startup, notices, confirmation, uncertain input and teardown.
- Existing `runtime_client.mjs`, `runtime_jobs.mjs` and
  `runtime_controls_browser.mjs` pass, including the new terminal/scientific
  boundary. Panel diagnostics: `/tmp/lcm-runtime-controls.pT9xwW/`.
- `runtime/test/run-bonito.sh` passes with the component in actual proxied Bonito
  gallery/deck-style frames and workbench hosts. Both themes, X-ray, binary Bonito
  sockets, independent runs, reconnect, surviving host and clean restart retain
  their prior gates. Latest evidence: `/tmp/lcm-terminal-relay.9ZPDrXWK/bonito-2.log`
  and `/tmp/lcm-runtime-bonito.1KK7yvEt/`.
- `runtime/test/run-broker.sh terminal`, `browser-tls-3.log`: **428 Julia assertions
  plus 15 real-browser checks**, exit 0. Browser keyboard/paste input crosses the
  real same-origin gateway, coordinator, TLS broker, agent and finite native Julia
  REPL. Verified variable evaluation, Unicode, multiline paste, session history,
  Tab completion, Ctrl-C, actual PTY dimensions, retained-state reconnect, clean
  namespace restart, 2,000 output lines, unaffected public health and no scientific
  jobs. Live screenshot: `/tmp/lcm-terminal-browser.q0dyhF/terminal-live.png`.
  The two-owner/private lease tests remain in the same suite.
- `browser-tls-1.log` also passed (initial live checks); its asynchronous child
  initially inherited null output streams. The test now forwards its bounded
  results to the harness. `browser-tls-2.log` records a failed multiline fixture
  that injected IME text rather than a browser paste event; the corrected test
  exercises xterm's actual paste handler, including Julia bracketed paste. No
  production deadline or authority rule was relaxed to obtain the pass.
- `component-2.log`: **101 component/runtime and 454 visual assertions**, exit 0.
  `publisher-1.log`: the complete publisher suite passes **1,133 assertions**,
  including presentation, ribbon, forms, uploads, workbench and X-ray contracts.
  The first ad-hoc invocation omitted the protocol import required by the test
  harness; its invocation was corrected, not the component logic.
- `publish-1.log`: `lcm playground build` completed all 21 Quarto documents;
  the terminal iframe and copyable source are present in `_site/widgets/index.html`.
  Owned JS syntax and `git diff --check` pass. The TLS harness removed its test
  broker and credentials; no `lcm.runtime.test=control-v2` broker remains running.
  Test processes were isolated; existing operator services were not altered.

**The complete P0–P7 goal remains active. This is verified progress, not release
completion.** Next: compose the same terminal into an actual registered Reveal
deck (the current test uses real deck-style Bonito frames), then implement the
shared scientific views and `CableStudy`/showcase consumers. Keep the verified
terminal and resource/lease/relay owners; do not rebuild them. Complete P6 and the
aggregate P7 gates. Actual quota-enforced Podman terminal launch and an actual
Docker runner remain external verification limitations on this machine, not
permission to weaken isolation. The host lacks delegated CPU control.

## P6 — registered scientific consumers and real browser composition

Implemented and verified on 2026-09-08:

- Added passive `StudyCases.LineParameters` and `CorridorImpedance` declarations,
  complete explicit specimen inputs, representative preparation inputs and
  bounded result-shape/unit projections. No solver, broker or UI dependency is
  imported by this case module. Authoritative worker validation remains unchanged.
- Added one shared `ScientificView(session, case, client)`: existing typed fields,
  `ScientificJob`, a persistent SVG plot, existing `DataTable`, worker selection
  and explicit preparation. Defaults derive from the case inputs rather than
  being repeated in the view. Incomplete/out-of-range edits disable submission;
  the existing job component owns provenance and stale-completion fencing.
- Added shared `CableGeometry` (local radial proportions only) and `StudyRuntime`
  (the existing selection/preparation/diagnostics components). One new structural
  stylesheet consumes the existing brand tokens in both themes. X-ray metadata
  belongs to the real components, not gallery replicas.
- Implemented the concrete `CableStudy` workbench and `Showcase` frame factories.
  Both use the same views, fields, job controls, diagnostics and Julia terminal.
  Both registrations consume one passive role declaration. Their fixed UI drivers
  are installed; no browser-provided command or implicit allocation was added.
- Replaced the showcase skeleton with a curated adaptation of all six legacy
  narrative sections, plus preparation and terminal demonstrations. Restored the
  original KU Leuven / EnergyVille SVG artwork from `5bc9916a`. Explicitly label
  the OHL/UGC curves as passive-length sensitivity, not statistical confidence
  bounds. Retain the original operation's solved/linearized reference assumptions.
- Added the shortcode's `requires-run=true` contract with a required public URL.
  Public slides immediately display their actionable placeholder without fetching
  an owned UI route. The actual registered Reveal deck now mounts the same terminal
  as the gallery and workbench; it is not merely a deck-style fixture.
- A shared `widget_shell(...; header=false)` option avoids duplicate headings for
  self-contained component frames. Scientific Run controls sit beside the inputs,
  not below the entire plot. A browser assertion checks their initial visibility.
  CableStudy's toolbar uses actual typed view-navigation actions.
- Added `SCIENTIFIC_CONSUMERS.md` and updated ARCHITECTURE. The small visible input
  set does not currently need file import/export; no persistent project service
  or implied REPL/model-cache sharing was introduced.

Issues exposed and corrected during verification:

- Corrected a workbench keyword-construction error in the first component run.
- Native numeric stepping is anchored at `min`: the earth-resistivity default
  was invalid with `min=0.01, step=1`. The corrected step and all default numeric
  fields now have explicit unit/browser regression checks.
- Screenshot review exposed HTML Observable wrappers escaping SVG text, and
  browser geometry checks exposed generic HTML property assignment failing to
  update SVG animated attributes. The shared views now use persistent native
  text nodes with `textContent`, and explicit `setAttribute` for circles/curves.
  No global Bonito monkey patch or vendored-package change was made.
- Corrected a fixture string comparison (`40` versus native `40.0`) to compare
  numeric values. The browser fixture now waits for actual Bonito socket/action
  bindings before its single navigation/input action, not merely DOM presence.
  No action retry or relaxed ownership/deadline was used to hide these failures.

Evidence:

- Full publisher suite: **1,239 assertions**, exit 0, latest
  `/tmp/lcm-scientific-consumers.3HONKlzE/publisher-layout.log`. This includes
  **96 scientific component assertions**, **462 visual assertions**, existing
  presentations, ribbon, toolkit, uploads, workbench and X-ray checks.
- Both actual profile environments ran `runtime/test/study_case_validation.jl`:
  **7 normalization/hash assertions per profile**, plus the dependency-boundary
  assertion. The complete UI inputs equal the workers' normalized inputs and
  produce identical hashes. This is not a numerical-execution claim.
- Runtime catalogue tests: **23 assertions**, exit 0,
  `/tmp/lcm-scientific-consumers.3HONKlzE/catalogue.log`. Installed registrations
  and passive definitions retain distinct launch behavior.
- `LCM_RUNTIME_BROWSER_SUITE=scientific bash runtime/test/run-bonito.sh`:
  **28 actual registered-consumer browser checks**, exit 0. Latest harness log:
  `/tmp/lcm-scientific-consumers.3HONKlzE/browser-final.log`; screenshots/PDF:
  `/tmp/lcm-runtime-bonito.75cu9Qnq/`. Both actual registered UI drivers ran behind
  the same-origin gateway. Verified public no-run behavior, independent drafts,
  valid/invalid input recovery, native SVG namespace/radius updates, persistent
  nodes, repeated light/dark parity, X-ray redaction, actual deck terminal mount,
  common immutable run prefix, overview return, and static print placeholders.
  The PDF has **11 pages, 1152 × 648 pt**, with no live frame requests in print.
- Earlier browser logs `browser-2` through `browser-5` retain the diagnosed
  failures; `browser-6.log` passed 26 checks before the final compact-layout gates.
  `browser.log` records an invocation error: the existing shell runner was not
  executable and is invoked through `bash`, not modified globally with chmod.
- `build-2.log`: all 21 Quarto documents rebuilt successfully. JS/shell syntax
  checks and `git diff --check` pass. Final test host directories are empty after
  cleanup; no operator service, broker/container or current Chrome session was
  changed. All test processes have completed.

**The P0–P7 goal remains active. P6 is in progress, not complete.** These tests
verify the real consumers and their passive contracts, not the final real
calculation gate. Next: drive the existing protected scientific job path from
both registered consumers, compare actual line/OHL-UGC results, and verify
cold/warm preparation, cancellation/retry, retained provenance and worker loss.
Reuse the existing TLS/job/artifact/lease harnesses and operation adapters; do
not rebuild them. Then finish aggregate P7 acceptance and handoff. Public
`assigned_execution` / `private_terminal` and terminal acceptance eligibility
remain gated; no native terminal fallback or isolation bypass was enabled.
The missing delegated Podman CPU controller and an actual Docker runner remain
external verification limitations, not a reason to stop useful P6/P7 work.

## P6 — first real scientific transport and projection hardening

Continued on 2026-09-08; this checkpoint does not close the P6/P7 gates.

- Found a real API/client mismatch: `RuntimeJob` required `assigned_execution`,
  but private control inventory never supplied it; only mock HTTP fixtures did.
  `control_capabilities` now supplies public discovery and private inventory from
  the same configured owners. It is not readiness or isolation evidence. Closing
  science/jobs removes execution availability. Private terminal acceptance stays
  separate and no physical preflight was weakened.
- Added `run-broker.sh scientific`, using actual registered Showcase/CableStudy
  hosts, browser controls, protected HTTP, role-authenticated TLS NATS, private
  MinIO and actual numerical profiles. The finite native driver is explicitly
  test-only. Added a separate direct-engine matrix comparison script.
- The first cold publication lost authority during initial transport setup.
  Job-connection startup now provisions/verifies configured scientific streams
  before reporting online; publication still rechecks policy. No lease or
  acknowledgement deadline was increased, and no numerical import moved into
  the coordinator. Initial stream provisioning is asserted by the new harness.
- A revoked unpublished job previously stayed displayed as queued because the
  browser coupled lifecycle polling to unavailable broker results. `RuntimeJob`
  now reads durable receipts independently, then fetches terminal data. Tests
  cover visible revocation without a result and recovery of delayed successful
  data without resubmission.
- Real results exposed a scientific display Observable inferred with
  `series::Nothing`. It could not accept a populated series. `ScientificDisplay`
  gives both states one type. The regression now renders the real component
  before delivering a first result, then invalidates it while retaining data.
- A non-sensitive `data-result-projection` acknowledgement distinguishes
  accepted/rejected/undelivered Bonito view updates. No input, result or private
  identity was added to X-ray metadata.
- Corrected two fixture issues: keep both launched UI consumers connected
  instead of allowing the unopened workbench's real disconnect grace to expire;
  inspect `data-result-current` on the actual inner job section, not its wrapper.

Evidence to this point:

- Shared capability/ownership gateway regression: **87 assertions**, exit 0.
- `runtime_jobs.mjs` and `runtime_client.mjs`: pass, including the new independent
  receipt/result failure cases. Updated the explicit Bonito HTTP fixture to
  follow the real receipt endpoint shape.
- Existing `run-bonito.sh`: pass with real binary sockets, shared scientific
  result bindings, terminal mounts, themes, X-ray, reconnect, independent hosts,
  stop/restart and cleanup; `/tmp/lcm-runtime-bonito.7hxopVTP/`.
- Full local runtime suite: **2,155 assertions in 115 testsets**, exit 0;
  `/tmp/lcm-runtime-regression.kX9LIK1i.log`.
- Scientific component regression: **102 assertions**, exit 0.
- Real TLS runs `qDQGshff` and `B4N43CVr` exposed lifecycle/publication issues;
  `Hv8e3Bw9` and `vuYfQGWC` completed the actual line calculation and exposed the
  result-to-display type error. Each path is under `/tmp/lcm-runtime-control.`.
- `/tmp/lcm-runtime-control.vK3HC3yO/scientific-live-failure.png` was inspected:
  real curves, numeric samples and provenance render correctly after the type
  fix. Its suite still failed because the assertion read the wrong DOM node;
  the subsequent fixture corrects that selector. Do not count it as a full pass.
- Some cold TLS retries exposed an existing MbedTLS buffered-reader close race
  (`ssl_unsafe_read` on a non-readable context) in vendored NATS `tls.jl:65`.
  This remains a recorded hardening item, not suppressed or claimed fixed.

The full publisher and corrected live-science gates are being rerun. Retain their
actual results below before marking any remaining acceptance item complete.
Cancellation/retry, replacement recovery, the second concrete consumer's real
numerical comparisons and aggregate P7/deployment checks still need evidence.

### P6 continuation — cold control compilation and TLS cleanup

- Full publisher regression completed: **1,245 assertions in 26 testsets**, exit
  0; `/tmp/lcm-publisher-regression.9X7AXTMH.log`.
- The corrected live run `Zby9ZesT` exposed a first-grant acknowledgement timeout.
  Fresh-process measurements found first-use lease JSON encode/decode work above
  the existing 2-second acknowledgement limit. Passive protocol precompilation
  reduced measured combined codec latency from about 2.74 s to 0.32 s (warm
  encoding about 20 microseconds). All **198 protocol assertions** pass.
- This was not the complete cold-path fix: `xqRMgzCf` passed assignment and
  preparation but revoked its first job. Opt-in compilation tracing in
  `run-broker.sh scientific` then identified a 1.62 s first-use scientific
  orchestration compilation burst in `lzbSEvKM`. That trace run lost preparation
  authority; its owned browser was stopped and the harness completed cleanup.
  Neither run is an acceptance pass.
- `StartupCompilation.jl` now requests compilation of installed orchestration
  entry points before live service announcement. It uses the existing managed
  scientific partition, not a new resource registry. It never invokes driver
  hooks, grants dummy leases, connects to transport or prepares a model. A new
  inertness regression passes **11 assertions**; the focused accompanying lease
  regressions pass **55 assertions**. No acknowledgement/presence deadline changed.
- The live browser fixture now fails promptly when preparation loses its
  assignment instead of waiting the full preparation deadline for an impossible
  ready state. `LCM_RUNTIME_TRACE_COMPILE=1` records an optional compiler trace
  inside the test's private temporary directory; it is not enabled by default.
- Fixed the recorded vendored NATS TLS reader close race with an exact guard:
  only MbedTLS's specific closed-read error on an actually non-readable context
  is normal cleanup. Other transport errors still propagate and the consumer
  buffer always closes. **13 focused assertions** pass; patch provenance is
  updated in `vendor/NATS/LCM_PATCHES.md`. Real TLS regression remains to rerun.
- `git diff --check` and the broker runner's shell syntax check pass.

The post-startup-compilation scientific run is in progress. Record its outcome
before closing P6. The full P0–P7 goal remains active.

### P6 continuation — acknowledged field edits and complete transport regression

- Full local runtime regression: **2,179 assertions / 117 testsets**, exit 0;
  `/tmp/lcm-runtime-cold-regression.nbH1rMka.log`. Includes the new TLS reader
  and passive startup-compilation tests.
- `sqnmNaGs` passed actual line-parameter and OHL/UGC calculation in the deck,
  then the workbench's initial line calculation. Its rapid edit-and-rerun gate
  found Run could capture the previous value before the newest field round trip.
  The job succeeded but the existing stale-input fence correctly refused to
  display it as current. This is a real interaction race, not a numerical failure
  or justification to discard stale-result validation. The whole suite failed.
- Added a shared ScientificJob input acknowledgement fence. A marked owned input
  scope immediately disables Run on native input/change, then accepts only its
  newest ordered Bonito acknowledgement and canonical draft. Older/duplicate
  replies cannot restore authority to submit an old draft. Display-only controls
  are excluded. Nested scopes remain independent and teardown removes listeners.
  Both scientific consumers use it; no domain validation was copied into JS.
- JS client/job regressions pass, including refusing Run while a local edit is
  pending. Real Bonito browser regression **passes**, with rapid numeric edits,
  an injected stale acknowledgement, immediate Run refusal and latest-value
  submission; `/tmp/lcm-runtime-bonito.RlD8ObpG/`. The subsequent expanded browser
  check adds checkbox/dropdown ordering and duplicate acknowledgement coverage.
- Scientific component tests: **102 assertions** pass. Focused RuntimeControls:
  **77 assertions** pass after correcting a test invocation missing the protocol
  imports normally supplied by `test/runtests.jl`; no production fix for that
  invocation error was necessary.
- The first full TLS rerun failed its initial cold PONG deadline (`pRe4Zf3r`,
  `/tmp/lcm-runtime-tls-regression.5m6M4V6z.log`). Explicit passive precompilation
  of the new TLS copy loop and existing buffered receiver corrected that path.
  The complete rerun **passes 306 assertions / 9 testsets**, exit 0:
  `/tmp/lcm-runtime-tls-regression.k4liniyx.log`, harness directory
  `/tmp/lcm-runtime-control.ocZZlI73/`. Includes separate authenticated identities,
  scheduler renewals/expiry, release/replacement, private TLS artifacts, real
  finite-child preparation/cancellation, and durable no-repeat recovery. No
  closed-reader warning occurred. Deliberate child-termination diagnostics are
  retained; they are not TLS failures. Harness-owned services were cleaned up.
- Added a real corridor cancellation/reprepare/retry browser gate after both
  consumers' calculations. It requires observed execution and acknowledged
  cancellation, retained last-good data, a distinct replacement executor and a
  fresh explicit result. Its finite whole-test budget includes one additional
  preparation window; production preparation/job/lease limits are unchanged.

The current numerical rerun is `/tmp/lcm-runtime-science-fence.YILM8gmc.log`, owned
harness directory `/tmp/lcm-runtime-control.JpIqLs1a/`. It has passed the deck's
line calculation and is preparing the corridor case. Cancellation/recovery and
the full consumer comparison are still pending; do not count them as passed.

### P6 continuation — retained numerical evidence and a valid cancellation gate

- `JpIqLs1a` subsequently completed all four real scientific cases, including
  the workbench's rapid edit/acknowledgement/rerun check. It did **not** pass the
  complete suite: the warmed 200-point corridor calculation finished before a
  status poll could observe execution, invalidating the new cancellation test's
  timing assumption. This is not evidence of failed production cancellation.
- Replaced that assumption with actual cold-preparation cancellation after an
  explicit release/reassignment. The gate requires worker-confirmed cancellation,
  retained last-good data, a usable lease, explicit preparation of a distinct
  executor and a new successful calculation. Running-job cancellation remains
  covered by the separate real TLS finite-child suite; these are distinct gates.
  No artificial solver delay, relaxed lease limit or unsafe execution fallback.
- Each completed numerical case now writes `scientific-partial-results.json`.
  Failure diagnostics also retain those results. Only the fully successful
  browser scenario writes `scientific-results.json`; partial data is never a pass.
- Expanded real Bonito browser regression **passes** in
  `/tmp/lcm-runtime-bonito.3uUozBj1/`, including native number/checkbox/dropdown
  ordering, old/duplicate acknowledgements and teardown.
- Full publisher regression **passes 1,245 assertions / 26 testsets**, exit 0:
  `/tmp/lcm-publisher-input-fence.trZBhu1c.log`. Full Quarto publication build
  passes, exit 0: `/tmp/lcm-publisher-input-build.vkTydBCt.log`.

The revised scientific scenario completed in
`/tmp/lcm-runtime-control.IbuYMDon/`, log
`/tmp/lcm-runtime-preparation-recovery.g2c2wUiS.log`. See the results below;
this run alone did not pass its enclosing Julia suite.

### P6/P7 continuation — conformance and cold-assignment investigation

- `IbuYMDon` passed **43 real browser checks**, including all four numerical
  cases, rapid input edits, preparation cancellation/replacement and worker-loss
  containment. Julia passed 29 assertions but failed an exact nominal-endpoint
  comparison: both consumers returned `150.00000000000003`, not exactly `150.0`.
  Consumer-to-consumer equality remains exact; the nominal endpoint assertion
  now uses Julia's ordinary approximate comparison. No engine code changed.
- Independent direct-engine validation passed **39 assertions**, log
  `/tmp/lcm-runtime-engine-reference.cmE8y19P.log`. The scientific broker runner
  now requires this separate numerical-environment check after the complete
  browser scenario; it never loads the engine into the coordinator.
- Registration-derived RuntimeConformance covers all nine current runtime and
  scientific component families through actual Bonito rendering and X-ray
  declarations: **411 assertions pass**. Full publisher regression including
  it passes **1,656 assertions / 27 testsets**. Legacy pre-v1 controls retain
  their separate tests; this is not a claim of full visual browser coverage.
- The aggregate acceptance runner records bounded, individually logged unit,
  browser, transport, scientific and host gates. Concurrent aggregate runs on
  this checkout are refused using a read-only lock on the script itself.
- Initial aggregate unit gates passed protocol **198**, execution core **205**,
  publisher **1,656**, worker **89** (including real power-flow preparation),
  and all three runtime JS suites. Its runtime test exposed a fixture depending
  on the caller's umask. The fixture now explicitly sets/asserts its deliberately
  broad permissions; production permission checks were not relaxed. The full
  private-umask runtime rerun passes **2,180 assertions / 117 testsets**, log
  `/tmp/lcm-runtime-private-permissions.QdGaZSud.log`.
- The initial aggregate shells were edited while running, causing a Bash
  continuation parsing failure. Their individual logs remain useful but those
  aggregate reports are invalid. Do not edit an executing shell harness.
- Fresh unchanged sequential aggregate `/tmp/lcm-runtime-platform.vRWQKvHG/`
  exited 1: publication/input checks passed, but the first scientific assignment
  failed before preparation. Earlier `PVbnWjhU` and subsequent traced
  `zopdsEsX` fixtures reproduce it. Concurrent heavy tests are therefore not a
  sufficient explanation. All fixture-owned services were cleaned up.
- The trace shows cold reservation/database specializations running while
  the inventory lock is held. The worker becomes stale before grant validation,
  which correctly refuses authority and reconciles the reservation. Startup
  compilation now additionally requests all placement and lease-control entry
  points. Inertness and cold-browser verification of this change remain pending;
  no acknowledgement, presence or lease deadline has been increased.
- Host report `/tmp/lcm-runtime-platform.f6hPnleW/results.tsv` exits 2:
  native/Podman effective CPU isolation is unavailable and `docker` is a Podman
  shim. No physical-isolation pass is claimed.
- Architecture, configuration and verification documentation now describe the
  owned process boundaries, conformance, aggregate commands and private SQLite
  backup/restore procedure. The latter was smoke-tested only on disposable
  databases in `/tmp/lcm-backup-check.JNrveFMY/`; operator state was untouched.

P6 remains open until cold assignment and the complete real scientific runner
pass together. P7 and the full goal remain active.

### Cold control-path follow-up

- The expanded startup compiler requests pass the 11 inertness assertions but
  alone did not fix the real gate: `HAYSsubh` granted revision 1, then released
  revision 2 after the first acknowledgement deadline. The trace still showed
  new database row/container specializations after the first row was inserted.
- Replaced query-dependent NamedTuple materialization with an internal `SQLRow`
  snapshot and a stable `Vector{SQLRow}`. Existing typed run/worker/lease/job
  conversion, SQL schema, ownership checks and deadline policies are unchanged.
  The added regression checks empty/populated/NULL/mixed-column results,
  detached values after statement/connection close and error recovery.
- Focused SQL, worker registration/migration, transactional allocation,
  multi-process capacity, agent leases, coordinator expiry/reconciliation and
  inert-startup regressions pass **221 assertions / 18 testsets**. One earlier
  focused invocation lacked the main runner's SQLite import; it was corrected
  without changing production code or tests to hide that invocation error.
- The browser test now fails immediately with the observed lease state/revision
  if initial assignment is lost, instead of exhausting a generic 90-second wait.
  `A638Bd7N` still failed the first ACK deadline, now with that precise reason.
  Query-dependent SQL container compilation was absent, but first nonempty
  diagnostic/occupancy collections still specialized inside the ACK window.
- Stabilized diagnostic response vectors, occupancy dictionaries and lease-ID
  sets across empty/populated states, and added acknowledgement-dispatch startup
  compilation. Focused strict-configuration, event-redaction/type, agent-lease
  and inert-startup tests pass **111 assertions / 7 testsets**.
  `X5ZQ391G` passed cold assignment, continued renewal and the real deck line
  result. Its power-flow child then exceeded the separate 120-second bootstrap
  bound while importing numerical packages. Leases continued renewing; the
  whole scientific suite failed and all fixture-owned resources were cleaned up.
- Profile modules previously imported their numerical engine before announcing
  the command reader, contradicting the bootstrap/preparation distinction.
  Added an optional, idempotent `load_profile!` hook (the four required hooks are
  unchanged). Both real profiles now defer fixed imports to explicit work in
  the child, report `loading_environment` during preparation and enter newly
  loaded methods through the existing latest-world execution boundary. No time,
  permission, resource or readiness limit was increased.
- Both actual profile input-normalization checks pass, including new assertions
  that registration/validation imports neither engine. The complete core suite
  passes **222 assertions / 12 testsets**, log
  `/tmp/lcm-runtime-deferred-core.RgMLKJM9.log`. Its 16 new framed-child assertions
  cover inert bootstrap/inspection, loading-stage reporting, fresh-method calls,
  warm cache reuse and physical cleanup after a deliberate loading failure.
- The scientific browser now fails promptly on an actual preparation failure;
  recovery distinguishes the new request from the preceding canceled request.
  The untraced complete scientific runner **passes**, exit 0, in
  `/tmp/lcm-runtime-control.Ze6jr9gp/`: **43 browser checks, 32 protected-runtime
  assertions and 39 separate direct-engine comparisons**. Four calculated-view
  screenshots and the complete `scientific-results.json` are retained. All owned
  children, test containers and credentials were cleaned up. Deliberate
  cancellation emits a child SIGTERM trace; this is expected test evidence.
- Measured cold/warm preparation: deck line **24.607 / 4.091 s**, workbench line
  **24.078 / 4.240 s**, deck corridor **178.569 / 4.392 s**, workbench corridor
  **168.317 / 4.097 s**. Explicit corridor cancellation/replacement preparation
  took **370.812 s**, within the unchanged 600-second profile deadline. These are
  measurements, not readiness promises or effective CPU/memory quota evidence.
- P6 is complete. The fresh sequential P7 unit acceptance group is running in
  `/tmp/lcm-runtime-platform.JU5ReNnV/`. The first permission review timed out
  before execution; its one permitted retry started normally. Broader regression
  and external host/isolation gates remain open; the full goal stays active.

### P7 — current sequential regression

The unchanged aggregate `unit` group **passes**, exit 0, with no unavailable
prerequisites: `/tmp/lcm-runtime-platform.JU5ReNnV/results.tsv`.

| Gate | Assertions / testsets | Seconds |
| --- | --- | --- |
| Protocol | 198 / 6 | 20 |
| Execution core | 222 / 12 | 45 |
| Runtime | 2,192 / 118 | 494 |
| Publisher, including registration-derived conformance | 1,656 / 27 | 170 |
| Worker, including actual PowerImpedance preparation | 89 / 10 | 283 |
| Runtime client / jobs / terminal JavaScript | All three suites pass | 1 / 0 / 1 |

Total Julia coverage: **4,357 assertions / 173 testsets**. Deliberate child
termination diagnostics are retained in the relevant logs; all gates exit 0.
This includes the fixed SQL row snapshots, cold-control container types and
deferred profile loading, not just the previously passing scientific fixture.

The sequential `browser` group is running in
`/tmp/lcm-runtime-platform.6NqBlQPA/`; it is not yet a pass. Visual inspection of
the successful scientific workbench screenshot found an unstyled footer link
using the browser's visited-link purple. The shared control contract now gives
plain links the existing `--lc-link` token at zero specificity, leaving owned
navigation/button/document styles authoritative. The scientific browser checks
this link through four alternating theme changes. No local palette was added.

Verification documentation now explicitly separates legacy Compose's optional
CPU overlay from mandatory v1 executor limits. The configuration handoff uses
the fully passing scientific run's timings rather than an earlier partial run.

The first current browser aggregate has exposed three failures; its independent
gates continue, and it cannot be counted as a complete pass:

- Gateway lifecycle: 86 control assertions passed but the first configured
  capability request exceeded its 10-second client deadline. Gateway startup now
  requests compilation of its exact configured HTTP entry before listening;
  the CLI also does so before starting broker scheduling. No route is invoked
  by this compiler request. Added three inertness assertions. Verification of
  this correction remains pending.
- Offline legacy publisher: `/widgets/form-toolkit` did not respond within its
  existing probe bound. The harness now copies child diagnostics to the retained
  acceptance log before cleaning its temporary directory. Cause and correction
  remain under investigation; no response deadline has been extended.
- Catalogue browser: the old P2 test expected five entries and an uninstalled
  scientific UI. P6 supplies six entries and both scientific consumers. The test
  now compares the public static catalogue with actual runtime registrations,
  requires showcase/CableStudy, and expects explicit scientific UI launch to be
  available independently of worker availability. Rerun pending. Its teardown
  test correctly found no runs because the obsolete assertion failed first.

Real Bonito browser integration **passes** in
`/tmp/lcm-runtime-bonito.2tcFRIxV/`. The registered scientific UI gate **passes 32
checks** in `/tmp/lcm-runtime-bonito.uSLndtvN/`, including the new link-colour
checks. Its PDF has 11 pages at 1152 × 648 points. Neither claims container quotas.

A separate test-harness audit found that preflight resolved a local engine but
subsequent shell calls could inherit remote connection settings. Added a
test-only launcher generator that reuses production detection, its exact local
command prefix and its filtered environment. Eight argument/connection-isolation
checks pass, including hostile inherited routing settings and quoted arguments.
An initial script-entry macro parse error was corrected. Aggregate integration
of the helper waits until the currently executing shell finishes; do not edit an
executing Bash harness.

### P7 — route readiness and cold-render investigation

The first browser group completed with **five failures**, not three:
`/tmp/lcm-runtime-platform.6NqBlQPA/results.tsv`. Presentation page loading and
ribbon theme initialization also timed out; the independent X-ray gate passed.
The second group, `/tmp/lcm-runtime-platform.OrM5j4wG/results.tsv`, also failed
five gates. Build, gateway lifecycle (**254 assertions / 11 testsets**), ribbon
and X-ray passed. Catalogue startup, two owned Bonito-host startups, the
scientific UI startup, presentation startup and the offline form GET did not.
Earlier complete scientific/P6 evidence remains valid; these current browser
regressions prevent a P7 completion claim.

Two actual readiness defects were distinguished:

- The offline probe previously used Downloads' default HEAD and repeated
  abandoned requests. It now consumes one real GET within the same 15-second
  bound and preserves failure diagnostics.
- The standalone publisher opened its listener before registering all routes.
  It now assembles the complete Bonito Routes value before constructing the
  server; a registration failure also closes an already-created broker client.
  Ribbon/X-ray fixtures follow the same ordering. Registration-derived unit
  checks were added. The second group's directory GETs all returned 200, rather
  than the intermediate 404 observed in the focused first-GET probe.

The form GET still exceeded its bound. A separate localhost-only compiler trace
in `/tmp/lcm-ui-cold-probe.log` measured package loading at **1.668 s**, first
form response **31.146 s**, and second response **0.052 s**. Trace output is
`/tmp/lcm-ui-cold-trace.jl`: the main costs are form/control rendering and
Bonito's initial page path, not a worker operation. A diagnostic -O1 run still
needed **24.309 s** for the first form; production compiler flags were not
changed. That run also followed a source edit and is not a controlled cache
comparison.

Disabled X-ray incorrectly called inspection hooks before checking policy.
The shared instrument entry now checks permission first, with a six-assertion
hook-invocation regression added. A network-free NoConnection/NoServer form
probe succeeded. A representative build-time UI workload using the already
installed PrecompileTools 1.3.4 is being evaluated; offline dependency resolution
changed no package versions. Its effect on actual HTTP startup remains unverified.

The aggregate engine wrapper is now integrated between runs. Its actual local
Podman transport group remains pending; the eight pure launcher checks are not
an engine-runtime certification.

Added operator-only `runtime/Caddyfile.example`, `proxy.example.toml` and
`PROXY.md` because the approved private proxy handoff lacked concrete templates.
They describe hashed two-user authentication, replacement/removal of asserted
identity headers, loopback forwarding and deployment denial checks. New tests
parse the real TOML with a disposable private key. **No Caddy installation,
configuration validation, service deployment or public exposure has occurred.**
Actual Caddy and effective Docker/Podman quota/remote-host gates remain open.

The UI workload package build succeeded in **59 s**. Repeating the original
HTTP probe with the same default compiler settings and trace options reduced the
first form response from **31.146 s to 8.458 s**; the second was **0.009 s**.
Evidence: `/tmp/lcm-ui-precompile.log`,
`/tmp/lcm-ui-precompiled-probe.log` and
`/tmp/lcm-ui-precompiled-trace.jl`. No production compiler flags changed.

The focused publisher architecture and X-ray regressions **pass 165 assertions**
(140 architecture, 6 permission-hook, 19 preview checks) in
`/tmp/lcm-ui-policy-tests.log`. The full offline publisher/graceful SIGINT gate
now **passes**, exit 0: `/tmp/lcm-shutdown-precompiled.log`. Every published
directory and every existing live probe returned 200 within its unchanged bound.
The gallery form took **3.275 s**, overlay toolkit **7.491 s**, data toolkit
**8.271 s** and repeater **5.628 s**. This verifies the latency correction on the
actual standalone publisher; registered multi-host browser reruns are ongoing.

The additional disconnected-workload test **passes 5 assertions** in
`/tmp/lcm-ui-inertness.log`. The actual proxy runtime configuration template
**passes 10 assertions**, including missing-key rejection, no storage creation,
operator-only administration and untrusted-peer denial. This is a TOML/runtime
identity test, not Caddy execution.

The focused catalogue browser gate now **passes**, with diagnostics retained in
`/tmp/lcm-runtime-catalogue.Fp96Yc2H/` and summary log
`/tmp/lcm-catalogue-precompiled.log`. Coverage includes 18 published routes in
both themes/desktop and compact flow, passive public/developer catalogues,
actual starter deck/template workbench/toolkit-gallery launches, preserved live
frame identity, native file upload/removal through the owned proxy and three
CLI teardown assertions. Worker availability does not gate static navigation.

A fresh aggregate browser group is running in
`/tmp/lcm-runtime-platform.DhDL45qP/`; do not treat its pending gates as passed.

The aggregate has now passed build (32 s), gateway lifecycle (160 s), offline
shutdown (65 s), catalogue browser (224 s) and Bonito browser (104 s).
Catalogue evidence: `/tmp/lcm-runtime-catalogue.2KlTcxkB/`; Bonito evidence:
`/tmp/lcm-runtime-bonito.5Jom44Yr/`. The latter includes independent hosts,
shared frames, actual binary WebSockets, reconnect, theme/X-ray identity,
explicit mocked scientific HTTP binding tests, owned host loss and clean restart.
Scientific UI, presentation, ribbon and X-ray gates are still pending.

Legacy broker/artifact harness cleanup now preserves the original failure and
makes owned-container cleanup failure fatal. Signal traps exit through that
single cleanup path. The artifact CLI helper has an exact test-owned name and
an ownership flag set only after successful create, so interrupted attach can
be cleaned up. The expanded fixture unit gate passes **18 assertions**:
8 argument/environment checks and 10 checks executing the actual cleanup bodies
against fake engine/file operations, including success, prior failure and
cleanup failure. No real container was touched by these unit checks. Actual
TLS transport and helper-container cleanup verification remain pending.

### P7 — complete current browser pass

The full aggregate browser group **passes**, exit 0, with no unavailable gates:
`/tmp/lcm-runtime-platform.DhDL45qP/results.tsv`.

| Gate | Seconds |
| --- | --- |
| Build | 32 |
| Gateway lifecycle | 160 |
| Offline publisher / graceful shutdown | 65 |
| Public/developer catalogue and owned UI launches | 224 |
| Real Bonito host isolation/reconnect/theme/X-ray | 104 |
| Registered scientific UI — 32 checks | 81 |
| Presentations, math notes, fragments, overview and PDFs | 178 |
| Ribbon/toolbar/interaction-state contract | 122 |
| Interactive X-ray, reset and redaction | 54 |

Scientific UI diagnostics: `/tmp/lcm-runtime-bonito.pfw3B4TK/`. The presentation
gate also checks all 18 published routes, both themes/compact navigation, the
actual math-note and native incremental-list interactions, and all three PDF
paths. This is the current combined regression pass, superseding the earlier
five-failure browser reports without erasing them.

The real transport group is running in
`/tmp/lcm-runtime-platform.tWGWIHpY/`. Its production-derived local Podman
launcher gate passed in 5 s. Runtime TLS, terminal TLS, legacy broker TLS and
legacy artifact TLS/cleanup are not yet complete. No effective executor quota,
genuine Docker or physical second-computer certification is implied.

The first fresh transport group completed with **two failures**:
`/tmp/lcm-runtime-platform.tWGWIHpY/results.tsv`. Terminal TLS **passes**
(155 s; 428 Julia assertions and 15 live Chrome/REPL assertions), retained in
`/tmp/lcm-runtime-control.vxaoUzpj/` and
`/tmp/lcm-terminal-browser.sqORBu/`. Legacy artifact TLS/role isolation and
exact named-helper cleanup **pass** (38 s). The pinned local Podman wrapper
therefore ran real container operations, not just argument tests.

- The legacy broker gate exposed a missing `import Bonito` in the newly extended
  architecture test. A focused command had imported it externally and masked
  the omission. The test now imports its own dependency; the standalone UI
  workload test does too. Full rerun pending.
- The first runtime control PONG missed its unchanged 2-second deadline.
  An isolated traced repeat reproduced this in
  `/tmp/lcm-runtime-control.kDi0jVKB/`. The configured NATS connection-controller
  closure spent **2.282 s** compiling before starting sender/receiver tasks.
  Sender compilation itself was only **0.053 s**. The vendored connection path
  now compiles that controller before returning the connection; it sends no
  synthetic request and changes no connection/lease deadline. Compilation
  failure closes its owned socket and completes the drain waiter.

Verification of that transport correction is running with optional compiler
timings in `/tmp/lcm-control-compiled-controller.log`. The broker harness's
existing trace option now applies consistently to full, terminal and scientific
groups. No production service, timeout, credential policy or retry allowance
was changed to obtain a pass.

The traced verification of the NATS controller correction **passes 306
assertions / 9 testsets**, exit 0, with diagnostics in
`/tmp/lcm-runtime-control.SJcNgb6Q/` and summary
`/tmp/lcm-control-compiled-controller.log`. The first direct BrokerControl
constructor now passes its original PONG check without fixture retries.
Coverage proceeds through actual authenticated role subjects, leases, targeted
durable delivery, production renewal/expiry, private TLS object storage,
scientific cancellation and result-before-ack recovery. SIGTERM traces are from
the deliberately canceled finite scientific children; no test failed.

A complete **untraced** transport rerun is active in
`/tmp/lcm-runtime-platform.ob68C9It/`; final unit regression follows it.
The browser pass above predates this small NATS scheduling correction, so the
final lifecycle/offline checks must also cover that current transport source.

## P7 — terminal admission order and retained startup diagnostics

The untraced transport aggregate `/tmp/lcm-runtime-platform.ob68C9It/` finished
with **two failures**, not a release pass: runtime TLS passed (306 assertions),
legacy artifact TLS passed, terminal browser readiness failed, and the legacy
worker missed its initial capability heartbeat. The missing test imports are
fixed and all publisher testsets before legacy transport passed.

The terminal failure screenshot/state is retained in
`/tmp/lcm-terminal-browser.r8cYtv/`. Its status displayed the initial selection
hint after losing usable assignment authority. The shared renderer now retains
the transport's explicit loss reason once a connection has been attempted.
The isolated Chrome/xterm regression passes **42 assertions**, including input
and caret disabling, visible assignment-loss reason and no automatic reconnect
or input replay. Evidence: `/tmp/lcm-terminal-status-browser.log` and
`/tmp/lcm-terminal-browser.nDIfwy/`.

An unchanged traced real terminal repeat passed all 15 Chrome checks and all
Julia testsets (`/tmp/lcm-terminal-traced.log`,
`/tmp/lcm-runtime-control.Tblqvg3R/`). This is not evidence that the intermittent
failure was fixed. The trace identified a fixture ordering hazard: it created a
second LocalIdentity gateway after admitting the browser lease, and that
gateway's first HTTP listener specialization took **3.001 s**, longer than the
unchanged **2 s** acknowledgement budget. Its read-only HTTP initialization now
finishes before the admission callback grants browser authority. Verification
is running in `/tmp/lcm-terminal-admission-order.log`. The gateway, browser and
lease deadlines were not extended.

Browser failures now retain current runtime inventory/assignments and browser
exceptions alongside their screenshot. The Julia fixture separately records
owner-filtered inventory, assignments and control events before cleanup.
The legacy worker failure was reproduced with an empty pre-connection log
(`/tmp/lcm-legacy-tls-diagnostic.log`,
`/tmp/lcm-legacy-worker-failure-5trksu/`). Its fixture now preserves owned worker
logs/process state and supports opt-in compiler timing using the same worker
module/arguments; ordinary acceptance continues to use the actual `lcm` CLI.
Root cause and a successful fresh aggregate remain required.

The gateway-before-admission repeat correctly rejected stale worker presence
(`/tmp/lcm-terminal-admission-order.log`); no stale reservation was granted.
Fixture admission now waits, within the existing 15-second setup bound, for a
new genuinely online challenged report. The next full terminal run **passes**:
103 protected gateway/lease assertions and 15 real Chrome checks, plus all
preceding terminal testsets. Evidence: `/tmp/lcm-terminal-fresh-admission.log`,
`/tmp/lcm-runtime-control.9cExTQVF/`, and
`/tmp/lcm-terminal-browser.PbIyda/`. A fresh aggregate repeat is still required.

Legacy compiler diagnostics reproduced its silent startup miss while the exact
owned worker was still running. `/tmp/lcm-legacy-worker-failure-kFwzya/` records
23.405 seconds of compilation, including 7.689 seconds for the CLI main entry;
the heartbeat bound remains 30 seconds. `LineCableModelsWorker` now passively
precompiles that entry while building its package cache. New import tests assert
no NATS connection and no scientific engine import. Compilation/legacy lifecycle
verification is pending; this is not yet a claim that the regression is fixed.

The current aggregate unit regression is active in
`/tmp/lcm-runtime-platform.uT8PojIH/`. No user service/browser, existing container,
operator configuration or isolation requirement has been changed.

The fresh aggregate unit run **passes every selected gate**, exit 0:
`/tmp/lcm-runtime-platform.uT8PojIH/results.tsv`.

| Gate | Assertions / testsets | Time |
|---|---:|---:|
| Fixture engine and cleanup | 18 / 2 | 9 s |
| Protocol | 198 / 6 | 23 s |
| Execution core | 222 / 12 | 49 s |
| Runtime | 2,202 / 119 | 539 s |
| Publisher | 1,699 / 29 | 149 s |
| Worker, including real PowerImpedance | 92 / 11 | 328 s |
| Runtime client / jobs / terminal JavaScript | All three pass | 1 s each |

Total: **4,431 Julia assertions / 179 testsets**, with no failed or unavailable
selected gates. The new worker inert-import checks pass. Deliberate SIGTERM
traces come from cancellation/recovery fixtures, not failed assertions.

A fresh complete transport aggregate is now running in
`/tmp/lcm-runtime-platform.jrOHKlDJ/`. The successful unit gate does not substitute
for the pending legacy CLI startup timing, combined transport or external host
verification gates.

## P7 — local regression complete; host-level release gates remain

The complete current transport aggregate **passes**, exit 0, with no failed or
unavailable selected transport gates:
`/tmp/lcm-runtime-platform.jrOHKlDJ/results.tsv`.

| Gate | Evidence | Time |
|---|---|---:|
| Pinned local fixture engine | Pass | 6 s |
| Runtime TLS | 306 assertions / 9 testsets | 277 s |
| Terminal TLS | 439 assertions / 16 testsets, plus 15 real Chrome checks | 165 s |
| Legacy broker TLS | 1,699 publisher assertions plus 49 lifecycle/authorization/mTLS assertions | 295 s |
| Legacy artifact TLS | S3 round trip and role isolation pass | 33 s |

Current terminal artifacts: `/tmp/lcm-runtime-control.xzPhe4NO/` and
`/tmp/lcm-terminal-browser.MZi25y/`. Current runtime transport diagnostics:
`/tmp/lcm-runtime-control.gKeS1Fgv/`. Owned test brokers/artifact containers and
temporary credentials were cleaned up by their passing harnesses.

The real untraced `lcm worker start` path now passes its previously failing
30-second startup heartbeat, replacement-worker admission, broker loss/reconnect,
queued cancellation, caching, expiry, child cleanup and graceful shutdown.
Neither that deadline nor the terminal's lease/acknowledgement bounds changed.
The corrected terminal admission order passes twice consecutively, including
the complete aggregate; the earlier intermittent failure is retained above.

Final current-source follow-up checks also **pass**, combined command exit 0:

- `/tmp/lcm-final-gateway-lifecycle.log`: 254 assertions / 11 testsets, including
  owner isolation, child/coordinator death, recovery, protected API and stalled
  broker startup/health behavior.
- `/tmp/lcm-final-offline-shutdown.log`: every published index alias and all 16
  widget probes return HTTP 200 with NATS absent, followed by graceful SIGINT.
  Cold form render is 3.152 s; the slowest toolkit probe is 7.540 s, below the
  unchanged 15-second request bound.

P3 is complete on this evidence. P4/P5 physical-resource acceptance and P7
overall release certification are **not complete**. The final read-only host
report `/tmp/lcm-runtime-platform.cE36Iwgu/results.tsv` exits **2** and records:

- Native user-systemd: **UNAVAILABLE — cpu_controller_missing**.
- Rootless Podman: **UNAVAILABLE — cpu_controller_missing**.
- Docker: **UNAVAILABLE — executable is a Podman shim, not Docker Engine**.
- Caddy is not available on PATH; its supplied operator template has not been
  rehearsed through an actual authenticated proxy deployment.

No host resource was started by these prerequisite checks. An approved host or
runner with effective CPU/memory/PID controls, a genuine Docker target, and the
approved proxy/remote-host context are required for the remaining gates. Do not
enable an unconfined terminal, weaken the resource contract, install Docker here,
or call finite native fixtures physical-isolation evidence. The next step is
operator direction for those targets, not another repeat of already-green local
suites. The goal remains unfinished pending that external access/configuration.

### Blocked checkpoint — 2026-09-08

The same external-access condition has persisted through three consecutive goal
turns. The latest read-only check still finds only `memory pids` delegated to the
user service, the Docker-to-Podman shim, and no Caddy on PATH. Existing local
alternatives, including create-only stopped-container policy checks, have already
been exercised; they cannot prove running resource isolation or another host.
The goal is now **blocked, not complete**, awaiting an approved host/runner and
access/deployment context. No test or service remains running for this audit.
Resume the remaining P4/P5/P7 acceptance work when those prerequisites change;
do not rerun green local suites solely to manufacture progress.

### Codespaces / Docker feasibility — 2026-09-08

User-authorized scope: a disposable proof of concept; **full Docker integration
remains deferred**. See [the reproducible report](runtime/CODESPACES_FEASIBILITY.md)
and [measured evidence](runtime/test/codespaces_feasibility_2026-09-08.json).

A fresh 2-core Codespace supplied genuine Docker/Moby 29.7.2-2 with cgroup v2.
The existing isolation verifier passed inside the pinned Julia 1.12.7 image for
plain and TTY-attached execution: 0.5 CPU, 512 MiB RAM, zero extra swap, 64 tasks,
8 MiB scratch, unprivileged identity and the required namespace/mount/security
checks. CPU throttling was observed. The missing-CPU negative control failed
closed with exit 78, Docker-level termination worked, and cleanup restored the
empty container inventory.

The environment has no real user-systemd manager, so current managed-agent
supervision/recovery is unavailable there. The initial probe also exposed the
production recipe/policy's `/usr/local/bin/julia` assumption: the pinned base
image's executable is `/usr/local/julia/bin/julia`. Record this for the deferred
integration; only the probe was adjusted. A probe stdin-EOF issue was isolated
and corrected without changing the application's terminal driver.

No application checkout, broker, gateway or worker service was deployed. No
production isolation policy changed. Scientific image builds, browser REPL,
resource-exhaustion/flood tests and managed-agent/lease/proxy/remote-topology
acceptance remain unproven by this pass. The obsolete Codespace was deleted with
explicit authorization; the disposable test Codespace was removed after its
results were retrieved. This completes the narrow feasibility request, not the
original full-platform goal.

### Both-engine acceptance resumed on owned Kubuntu — 2026-09-08

The user superseded the Docker deferral and authorized scoped installation,
sudo-backed host configuration and end-to-end Docker plus Podman verification
on `ts ssh kubuntu`. This alias is handled by the user's rootless Tailscale CLI;
it invokes OpenSSH through `tailscale nc`, not a direct ordinary SSH hostname.

The target identifies as `amauri-a70mob`, Ubuntu/Kubuntu 24.04.4, kernel
6.8.0-139-generic, account `amauri` / UID 1000. Its real user-systemd manager
already delegates `cpu memory pids`; subordinate UID/GID ranges are present.
Pre-existing failed desktop units were inspected and left untouched. No
container engine, engine configuration or container data existed before setup.

The user entered sudo's hidden password prompt directly in a separate local
Konsole running `ts ssh -tt kubuntu`; no password traversed chat or was captured
in logs. The reviewed root bootstrap installed official Docker CE 29.8.0,
rootless extras/Buildx/Compose and Ubuntu Podman 4.9.3 with crun, uidmap,
slirp4netns and fuse-overlayfs. Newly installed rootful Docker/containerd services
were stopped and disabled. Docker now runs as a user service in rootless mode;
Podman is rootless too. User lingering was enabled. No docker-group membership,
passwordless sudo, firewall weakening or local Docker installation was added.

Both engines report cgroup v2 and CPU/memory/PID controls. Their initial container
inventories were empty. The shared container policy, recipe and test fixture now
use the pinned base image's actual `/usr/local/julia/bin/julia`. Focused local
tests pass **260 assertions / 4 testsets**, including recipe/launch parity for
both engines. Image builds cap Julia precompile tasks at two via a build-only
ARG, without adding an undeclared runtime environment variable.

Private acceptance root:
`/home/amauri/lcm-runtime-acceptance/run.Rp0IxHde/`. The allowlisted current-source
snapshot excludes ignored state, credentials and Git metadata; archive SHA-256
`13104b2048327fe7ccc7f78bce7639dd7ef324b9579952c68eebcd1be0f2818e`.
The later build-only ARG was synced separately. Julia 1.12.7 was downloaded from
the official distribution and verified against SHA-256
`4e7e9e776634d24835250de67cde39b0d4af15bc432eb20697e6be6c28ea69e8`.
It and its depot are confined to the private acceptance root, without changing
the remote user's existing shell or toolchain. Runtime preparation and the first
real Docker terminal-image build are underway; **no physical/end-to-end pass is
claimed yet**. Build and setup logs are under the acceptance root's `logs/` and
`/tmp/lcm-host-provision.tjcWDntY/`; local provisioning artifacts are in
`/tmp/lcm-kubuntu-provision.pmG1bOOq/`.

### Physical terminal isolation passes on both engines; managed relay pending

The production journal, container policy, PTY driver and in-container kernel
guard now pass the opt-in `runtime/test/physical_terminal.jl` on the owned host:

| Engine | Assertions | Log under the acceptance root |
|---|---:|---|
| Rootless Docker | 23 / 23 | `logs/docker-physical-terminal-r4.log` |
| Rootless Podman | 23 / 23 | `logs/podman-physical-terminal-r5.log` |

These exercise two private REPLs, exact ready-marker admission, output flooding,
kernel `memory.events` OOM-kill evidence, process-limit exhaustion, surviving
unrelated session state, and exact receipt/container retirement. Inventories
return to baseline. This is real physical isolation evidence, not a native
fixture; it does **not** certify the managed gateway/browser path.

The immutable terminal references are Docker
`localhost/lcm-julia-terminal@sha256:b51bc0192278a53ad70dea81ee7047e150007ad99dc7227a8dd0cf1a7d14cda8`
and Podman
`localhost/lcm-julia-terminal@sha256:ef5cf0af4503c2fb83aa2817175b2d8b8b7deed2ef1e79ce1dbed2d9e524d596`.
Both are built from the same current shared guard and pinned Julia base.

Actual engine observations required narrow normalization: Podman's bare full
image ID and scalar single entrypoint, explicitly configured `HOSTNAME`, its
supplementary finite `RLIMIT_NPROC`, and rootless OCI default device ownership.
Only the six exact standard device numbers and three exact read-only metadata
files are admitted; arbitrary environment values, mounts and relaxed limits are
still rejected. Kernel guard unit tests pass 120 assertions. A full local runtime
regression passed before the final supplementary-NPROC/cleanup changes; focused
follow-up tests cover those latest changes separately.

The new `run-broker.sh physical` gate uses an actual verified user-systemd agent,
TLS NATS, two owner-fenced assignments and the private WebSocket gateway. Its
first attempts fail on the first terminal action with WebSocket closure 1006,
after 13 successful checks. Redacted receiver/watchdog diagnostics are being
used to locate the first failure. Real service stop also reached systemd's
60-second deadline; successful post-stop recovery is not counted as graceful
shutdown. Scheduler-interruption cleanup has a regression under verification.
All failed runs retained diagnostics and cleaned their exact owned resources.

The standalone Caddy 2.10.2 test binary was verified against the official SHA-512
manifest and installed only inside the acceptance root. No proxy listener or
public deployment was started. Scientific image builds and actual proxy,
managed terminal/science, remote-topology and browser acceptance remain pending.

### Actual Docker proxy/terminal and graceful service stop — verified

The untraced separated-client gate passes through the shipped Caddy template with
verified TLS: **126 private-client checks + 29 managed/proxy checks** in
`logs/docker-physical-proxy-separated-r3.log`, exit 0. Diagnostics are under
`/tmp/lcm-runtime-control.Pqh4ZLeX/physical-agent/`. Anonymous/forged identity,
direct private-route bypass, researcher-to-administrator spoofing and foreign
Origin checks reject access. Two distinct owner REPLs, input, restart, private
state, keepalives, retained state, graceful stop, socket revocation and exact
container/journal cleanup pass. Caddy's listener and certificate trust are private
to this rehearsal; no public hostname or global CA trust was configured.

An earlier unproxied Docker run also passed **147 checks**, including deliberate
SIGKILL and post-stop recovery with two live containers, in
`logs/docker-physical-agent-trace.log`. Earlier failed runs remain recorded above.

Two corrections matter for interpreting those failures:

- The initial combined-process test let client TLS/JIT compilation pause its
  coordinator's short lease/ACK windows. The physical harness now launches
  `physical_terminal_client.jl` separately, matching the actual browser/process
  boundary. Both clients send serialized keepalives during startup and long
  numerical work. Terminal codec/request specializations compile passively before
  advertising control availability. No lease, ACK or socket deadline was extended.
- Julia's SIGINT-as-exception could land in an unrelated scheduler; trying an
  exit hook was not sufficient for joined asynchronous cleanup either, and that
  experiment was removed. The root now also observes its already-verified
  systemd incarnation for `deactivating` and closes cooperatively before process
  exit. A replaced/missing unit fails closed; executor admission still requires
  `active`. Failed scheduling tasks cannot skip physical retirement. This path
  passed actual graceful stop, not only unit tests. Exact typed supervisor tests
  pass **58 assertions**; the latest focused policy/cleanup/layout suite passes
  **403 assertions**. A fresh complete runtime regression is running separately.

Final Docker numerical images built successfully: line-parameters
`localhost/lcm-line-parameters@sha256:4f4d8a9667a090c2669c0d1b85391c1bfe12d94b093c9a0a9f85471b88fa9968`
and power-flow
`localhost/lcm-power-flow@sha256:ded77f3da6be19e04697586c91ac46acedfce89df4a4ad55efc959b08f439e32`.
Podman line-parameters is
`localhost/lcm-line-parameters@sha256:062c9f4c3feae4014955886640b682664b007b9b78312ee4a4b0e758e5f0ec7c`.
Build success alone is not numerical acceptance. The new opt-in
`physical_science.jl` checks both real scientific profiles through protected
preparation/job APIs while the private terminal client remains active. Its runs,
Podman managed/proxy acceptance, physical cross-computer topology and final
browser/full-release certification remain pending.

### Fresh pinned-engine inspection and repeatable physical acceptance

The complete local runtime regression finished successfully: **2,256 assertions /
121 testsets**, exit 0, `/tmp/lcm-runtime-remote-host-regression.log`.

The first Podman managed line/terminal run
`logs/podman-physical-proxy-line-r2.log` failed its readiness gate (1,226 passing,
two failed assertions and one consequent missing-receipt error). It retained two
live REPLs and a prepared numerical process but could not obtain fresh ready
evidence. Reading its protected status returned `idle/unknown`, with each
inspection taking approximately five seconds. A direct `podman info` measured
0.67 seconds: live pre/post checks repeated engine discovery and four info calls
each, exceeding the existing five-second report window. The test cleaned its
exact resources and service; stale readiness was not accepted as success.

`recheck_container_host` now reads the already pinned local command once and
derives both fresh prerequisites and exact engine/storage identity from that
response. There is no prerequisite cache, implicit context switch, relaxed
deadline, or skipped per-container/lease/kernel policy. New negative cases cover
replaced storage/daemon identity, missing controllers, malformed/failed responses
and remote command rejection. Host, recovery and policy tests pass **340
assertions / 11 testsets** in `/tmp/lcm-runtime-host-recheck-r2.log`. The first
attempt at this new test had a fixture declaration-order error; it was corrected
before the passing run. A complete regression is also rerunning after this change.

The corrected Podman managed gate is running in
`logs/podman-physical-proxy-line-r3.log`; the remaining Podman power-flow image is
building in `logs/podman-power-flow-build-r3.log`. Neither is yet counted as passed.

The aggregate now has an explicit `physical` group, requiring all three installed
manifest-pinned images and Caddy before allocation. It runs actual kernel/PTY
limits and managed science/two-owner terminal checks for both graceful and crash
stop modes, selecting Docker or Podman through the existing engine resolver.
It accepts an exported checkout and does not require Node for this non-browser
group. Syntax and missing-image fail-closed checks pass. Usage and the remaining
verification boundaries are documented in `runtime/PHYSICAL_ACCEPTANCE.md`.

### Both numerical images available; Docker full graceful gate passes

The corrected Podman line/REPL/proxy run passed, exit 0: **222 separate-client +
109 managed checks**, `logs/podman-physical-proxy-line-r3.log`, diagnostics
`/tmp/lcm-runtime-control.2cdPO0oC/physical-agent`. Fresh inspection took 1.753
seconds; the real line result, retained REPL state, agent SIGKILL and exact
post-stop recovery all passed. The latest host-recheck full runtime regression
also passed **2,275 assertions / 122 testsets** in
`/tmp/lcm-runtime-remote-host-regression-r2.log`.

Podman power-flow finished building with the same guarded recipe:
`localhost/lcm-power-flow@sha256:9622426c590dfacccd7bdc550804f57b1b83e8030c62817e3383c4abd54cf76e`.
All three approved images now exist on both engines; build logs remain private
under the acceptance root. Build success is not an execution result.

The first full Docker aggregate, `/tmp/lcm-runtime-platform.q7SHtxTC`, passed
engine/host/kernel gates but failed both managed modes at the first scientific
job. Durable timestamps show a submitted job becoming revoked while its lease
was released. Concrete managed job-delivery functions were missing from the
passive startup compiler requests. Additionally, scientific HTTP traffic still
ran inside the test coordinator, unlike an actual browser. Both were corrected:
all test-client actions now run in the separate client process and acquire their
scientific assignments through the protected API; managed delivery/admission/
persistence paths compile before advertising availability. No synthetic job,
grant or model is run to warm them, and all control deadlines are unchanged.
The complete regression after that change passes **2,277 assertions / 122
testsets**, `/tmp/lcm-runtime-job-compilation-regression.log`.

The revised Docker aggregate is `/tmp/lcm-runtime-platform.wsc7p9pe`:
its full **graceful** managed gate passes **410 client + 30 managed checks**,
both real line-parameter and power-flow results, with fresh inspections of 0.339
seconds. Log: `physical-managed-graceful.log`; diagnostic root:
`/tmp/lcm-runtime-control.4APew4aH/physical-agent`. Its forced-crash mode is still
running; no aggregate pass is claimed yet. Fixed debug records now distinguish
missed acknowledgement/authority timing without logging commands or credentials.
Post-stop lease expiry in this successful test is expected, not new authority.

A separate `physical_consumers.jl` gate is staged for the actual local
Bonito/Reveal browser consumers against separately provisioned remote agents and
private S3 storage. It reuses `scientific_live_browser.mjs`, including an opt-in
large-result S3 retrieval check; that browser check alone makes no OS-isolation
claim. Syntax checks pass; the physical cross-computer run has not started.
Private staging directories are `/tmp/lcm-physical-consumers.N4E6lsnP` locally and
`/home/amauri/lcm-runtime-acceptance/consumers.RiWNuALa` remotely. Provisioning
helpers are test/operator actions, not a new SSH runtime orchestrator.

### Complete Docker physical matrix; Podman repeat in progress

`/tmp/lcm-runtime-platform.wsc7p9pe/results.tsv` exits **0**: all five gates
pass, with no unavailable prerequisite. The forced-crash mode passes **405
client + 28 managed checks**, including both actual scientific profiles and
independent REPLs through Caddy; its diagnostic root is
`/tmp/lcm-runtime-control.C3ONZW6g/physical-agent`. Graceful stop took 280 seconds
and forced-crash acceptance 278 seconds overall, including cold startup and
numerical preparation. These are full test durations, not preparation timings.
Both modes restore exact owned-resource inventory and retire lease authority.

The first Podman aggregate `/tmp/lcm-runtime-platform.SxYBSFih` passed its host
and 23-check kernel/PTY limit gates, then failed before either managed launch:
the pinned NATS fixture image was absent from Podman's separate image store.
It was explicitly pulled; no implicit image download or relaxed executor policy
was added. The repeat is `/tmp/lcm-runtime-platform.vwGQYzSE`, with broker and
executors all on rootless Podman. Its managed modes are still running.

The cross-computer consumer gate now also reuses the actual Chrome/xterm
terminal check after scientific worker loss. It keeps the existing renderer,
wire protocol and resource policy, including Unicode, multiline input, history,
completion, interrupt, resize, reconnect, restart and both themes. This is still
a staged gate, not a passed result.

### Both physical engine matrices pass

The full rootless Podman repeat `/tmp/lcm-runtime-platform.vwGQYzSE/results.tsv`
exits **0**, with all five gates passing and no unavailable prerequisite. Its
graceful managed mode passes **499 client + 30 managed checks** (316 seconds);
the forced-crash mode passes **499 client + 28 managed checks** (317 seconds).
Both prepare and execute the real line-parameter and power-flow profiles while
two private terminals remain independent. Both actual NATS fixture and executors
use Podman. The Caddy boundary is the same one used for Docker acceptance.

The completed Docker and Podman runs restore the whole-container baseline to
empty on each engine and unload the exact `lcm-agent-worker-a.service`; no stale
owned executor, lease authority or scratch remains. Scope guards, kernel limits,
the 2-second ACK and 5-second science freshness bounds are unchanged. Local
passive client regressions also pass: terminal transport 47 assertions and the
scientific job projection/fencing/cancellation suite.

The separate-computer browser/S3 rehearsal is now starting. The engine matrices
are complete; browser cross-host acceptance is not yet claimed.

The remote S3 fixture setup required two explicit corrections before any agent
or consumer run: the exported test checkout lacked `artifact_initialize.sh`
(the bind option created an empty directory, which was verified empty and
removed), and MinIO then hit the desktop account's existing inotify-instance
limit. A fresh standalone `inotify_init1` also returned **EMFILE (24)** with
`max_user_instances = 128`, confirming this was already a host-user resource
limit, not a scientific executor failure. The fixture now preflights its script
and runs MinIO as a separate non-root subordinate UID with only its two private
TLS files mapped to that UID. No global sysctl, file limit, executor constraint
or unrelated application was changed. Failed setup containers and enumerated
secrets were removed; numbered non-secret logs remain in the rehearsal root.

The non-root UID experiment still returned EMFILE: the parent user namespace's
watcher accounting remained exhausted. No successful workaround is claimed for
that attempt. The rehearsal instead uses the already verified private Caddy
binary to terminate HTTPS in front of a literal-loopback-only MinIO backend;
runtime S3 clients retain CA/hostname verification. The test initializer gained
an explicit, fixed loopback HTTP mode (default remains verified HTTPS), not an
arbitrary URL or a plaintext runtime storage configuration. The Caddy listener
is private, has no admin listener, and is owned/retired by the fixture.

### First physical browser run exposed a terminal reconnect race

Docker's first cross-computer science run passed **49 browser checks**, including
both actual calculations in both registered consumers, cold/warm preparation,
cancellation/replacement, line survival after power-worker stop and a **95,261-byte
private S3 result** (not inline data or a shared filesystem). Evidence is in
`/tmp/lcm-physical-consumers.N4E6lsnP/browser-r1`. The complete gate nevertheless
failed **24 passed / 2 failed**: a test checked its remote-stop marker before the
operator process finished, and the real terminal browser failed at reconnect.
No complete pass is claimed for that run.

Intentional disconnect previously closed the browser socket immediately; the
gateway still owned its attachment while awaiting the remote writer-release
reply. Reconnect could race that cleanup and be rejected despite a valid lease.
The shared browser client now serializes `disconnect`, disables input/reconnect
while awaiting its acknowledgement, and retains a visible disconnecting state.
The gateway avoids a second disconnect after confirmed release. A new claim may
wait at most five seconds for an **already closing** predecessor, rechecking
current owner/lease and never stealing an active attachment or freeing an
unreconciled slot. Uncertain input is still never replayed.

Delayed-ack browser transport tests pass **56 assertions**; full runtime tests
pass **2,286 assertions / 123 testsets**, exit 0, in
`/tmp/lcm-runtime-terminal-handoff-regression.log`. A focused real remote Docker
terminal run then passed **15 Chrome + 11 coordinator checks**, including
reconnect/restart, in `terminal-docker-r2.log`; its screenshot is retained at
`/tmp/lcm-terminal-browser.T8G46A`.

The physical consumer fixture now keeps the actual presentation connected during
its companion terminal test, runs the private operator stop in a separate Node
process, waits for that process's result, and performs numerical parity assertions
after closing live coordination. Test-only Julia compilation must not consume a
live ACK budget. The full corrected Docker repeat is running in
`browser-docker-r2`; the old agent incarnations were stopped and their exact
journals reconciled before the fresh `docker-r2` configurations were started.

### Corrected Docker cross-computer consumer gate passes

The complete Docker repeat exits **0**, with **49 scientific browser checks,
15 actual Chrome terminal checks and 29 coordinator/ownership checks** passing
in 597.5 seconds. Evidence: `browser-docker-r2.log` and `browser-docker-r2/`
under `/tmp/lcm-physical-consumers.N4E6lsnP`; terminal screenshots are retained
at `/tmp/lcm-terminal-browser.YWRUz2`. Both real applications prepare and run
line parameters and power flow; cancellation/retry, independent worker loss,
warm reuse and private S3 retrieval pass. Both exact Docker agent incarnations
were subsequently stopped and their acquisition journals reconciled. Inventory
then contained only the two owned broker/artifact fixture containers.

The separate direct-engine comparison initially rejected eight near-zero
admittance coefficients: cross-machine cancellation residuals around 1e-27
were compared using relative-only scalar `≈`, including against exact zero.
The test now retains Julia's default relative tolerance and adds an absolute
floor of **one Float64 epsilon times the reference frequency-slice scale**.
There is no unit-dependent fixed tolerance or scientific implementation change.
Eight positive/negative controls check roundoff, meaningful near-zero errors,
relative errors, NaN and zero-scale handling. All **39 engine comparisons +
8 tolerance controls** pass; original and corrected logs are retained as
`parity-docker-r2.log` and `parity-docker-r2-bounded.log`.

Fresh Podman worker services are started with verified fixed recovery. Their
full separate-computer browser repetition has not yet passed.

The final runtime regression, including the five-second unresolved-attachment
deadline and ownership retention, passes **2,289 assertions / 123 testsets**,
exit 0 (`/tmp/lcm-runtime-final-remote-regression.log`). Shared terminal client
tests pass **56 assertions** again. A fresh local Podman TLS terminal regression
passes **435 assertions / 16 testsets**, including 99 protected live-REPL gateway
checks, exit 0 (`/tmp/lcm-terminal-final-tls.log`, fixture diagnostics
`/tmp/lcm-runtime-control.iwMbo5wZ`). Its broker and test credentials are cleaned
by the harness. This is transport/ownership evidence, not an assertion that the
IT-managed local host gained the missing CPU controller.

### Podman science passes; browser startup deadline corrected

`browser-podman-r1` passes all **49 scientific browser checks**, including both
consumers, cancellation/replacement, worker loss and the same **95,261-byte S3
artifact**. The complete gate nevertheless exits 1 (**28 passed / 1 failed**):
the terminal browser's generic DOM wait expired after about 7.5 seconds while
the real terminal was still starting. Its lease remained active/usable, the
worker was online and there were no browser errors. The shipped terminal
startup bound is **120 seconds**; the fixture had not used it.

The browser fixture now receives that existing bound from
`TerminalSessionLimits()` for initial startup and explicit restart only. Normal
interactions/reconnect retain 7.5 seconds; explicit failed/exited/uncertain/
disconnected states fail immediately. The parent and companion fixture waits
are bounded to cover those two allowed startup periods. No production timeout,
resource limit or startup policy changed. A focused real Podman repeat reached
ready in **9.640 seconds**, demonstrating the old test cutoff was premature;
its remaining interaction checks and the full repetition are still pending.

The two complete engine matrices' twelve selected non-secret reports/logs are
also retained remotely in `run.Rp0IxHde/logs/physical-engine-evidence.tar`
(SHA-256 `17c08e67bbadb5820fa7c1f15c4530361d0d548b0728d1a46f20fac4a05f4590`),
independently of the host's temporary-directory retention policy.

The focused Podman terminal repeat exits **0**, passing all **15 actual Chrome
checks + 11 coordinator checks** in 72 seconds (`terminal-podman-r2.log`,
screenshots `/tmp/lcm-terminal-browser.7MU1cM`). The corrected shared browser
fixture separately passes **42 rendering/interaction checks**, exit 0
(`/tmp/lcm-terminal-final-browser-fixture.log`, screenshots
`/tmp/lcm-terminal-browser.0bwaFf`). Podman's first-run science results also pass
**39 direct-engine comparisons + 8 tolerance controls** (`parity-podman-r1.log`).
Both first-run Podman agents were stopped/reconciled; fresh `podman-r2` services
now run the complete scenario again. No combined full-run pass is claimed yet.

## Final mandatory acceptance — complete, 2026-09-09

The complete corrected Podman run exits **0**, passing **49 scientific browser
checks + 15 actual terminal browser checks + 29 coordinator/ownership checks**
in **717.2 seconds**. Guarded terminal startup took **10.173 seconds** within
the unchanged 120-second production startup bound. Its separate engine process
then passes **39 numerical comparisons + 8 tolerance controls**, exit 0.

Final local evidence root: `/tmp/lcm-physical-consumers.N4E6lsnP`.

| Engine | Full combined log / evidence directory | Direct-engine log |
|---|---|---|
| Rootless Docker 29.8.0 | `browser-docker-r2.log`, `browser-docker-r2/` | `parity-docker-r2-bounded.log` |
| Rootless Podman 4.9.3 / crun 1.14.1 | `browser-podman-r2.log`, `browser-podman-r2/` | `parity-podman-r2.log` |

Podman's final terminal screenshot is `/tmp/lcm-terminal-browser.UAfvR7`.
Both consumer runs use actual Bonito/Reveal on the gateway computer, separate
managed agents on owned Kubuntu, immutable executor images, verified TLS NATS
and private HTTPS S3. Both retrieve the actual **95,261-byte** result through
the owner-authorized artifact endpoint, while the public digest route denies
access. The cross-computer fixture shares Docker-hosted broker/storage services
between engine repetitions; the separate complete Podman physical matrix uses
Podman for its broker and executors. No shared result filesystem stands in for
the transport, and the coordinator never imports Bonito or scientific engines.

Measured preparation includes browser status polling, not an SLA:

| Engine | Line cold, deck / workbench | Power-flow cold, deck / workbench | Warm explicit preparation | Cold cancellation/retry recovery |
|---|---|---|---|---|
| Docker | 26.608 / 24.239 s | 108.429 / 108.655 s | 4.121–4.134 s | 110.680 s |
| Podman | 30.596 / 30.221 s | 120.559 / 130.641 s | 4.122–4.137 s | 132.806 s |

Scientific executors retain 1 CPU, 4 GiB memory, 256 PIDs and 128 MiB scratch;
the terminal retains 0.5 CPU, 512 MiB, 64 PIDs and 8 MiB scratch. Engine/kernel
acceptance, owner fencing and lifecycle deadlines were not weakened for a pass.
The latest runtime regression remains **2,289 / 123 testsets**; final terminal
TLS is **435 / 16**, shared transport **56**, and isolated Chrome rendering **42**.

Final read-only audits exit **0**:

- Both remote engines are rootless and their complete container inventories are
  empty. All eight rehearsal agent incarnations have stopped markers and exact
  journal reconciliation; neither worker unit exists. Only inert journal owner
  metadata/locks remain, with no resource receipt or executor scratch.
- The two owned fixture containers, initializer and private Caddy proxy are gone.
  TLS/HTTPS fixture listeners are closed. Enumerated remote credentials,
  certificates and agent configurations are removed; diagnostics are retained.
- The exact local SSH forwarding process and both forwarding ports are gone.
  All six local rehearsal UI-host scratch roots are empty, and enumerated local
  credentials, certificate files and private stop-command configurations are gone.
- Rootless Docker's user service remains active. Newly installed rootful Docker,
  its socket and containerd remain disabled/inactive. Podman, approved images,
  private tools/depot/checkout and non-secret acceptance evidence remain installed.
  No unrelated workload, global inotify limit, local engine installation or
  user browser profile was changed by this rehearsal.

All mandatory P0–P7 work is complete. The earlier failed attempts remain in this
chronological ledger; they are not counted as passing runs. Production public
hosting/accounts/certificates, project CRUD, HA, arbitrary hostile-code hosting
and optional workbench file import/export remain outside the locked v1 scope.

Final handoff integrity: all **140** runtime/protocol/common/execution-core/
worker/broker source files match between local and remote checkouts, verified
before synchronizing the final client/test/documentation updates. The source
manifest SHA-256 is
`5e96ce63b80cd6ba0ae4b2a1557b80ba91cb9931f48c2d2be288a10f77d5cb9a`.
Selected successful cross-host result files, screenshots, regression logs and
the final audit script are retained in the remote private acceptance root at
`run.Rp0IxHde/logs/cross-host-evidence.tar`, with the same local/remote SHA-256
`c1538e5ea1aa4f883923f4b3f63adf9194ff4aa58896268004321715b730e3f7`.
No credentials, complete browser profile or runtime database is included in that
archive. Final source syntax and `git diff --check` pass. No commit, branch
switch, registry push or public service deployment was performed.
