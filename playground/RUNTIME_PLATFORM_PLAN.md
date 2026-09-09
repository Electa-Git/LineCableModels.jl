# Runtime platform v1 — execution plan

Status: complete — all mandatory P0–P7 gates verified, including both rootless
engines and separate-computer consumers, 2026-09-09. Execution evidence is tracked in
[RUNTIME_PLATFORM_PROGRESS.md](RUNTIME_PLATFORM_PROGRESS.md).

Prepared on 2026-09-07 against playground commit `d50367bf`.
Legacy content reference: `feat/codespaces-showcase` at `5bc9916a`.
Recheck the worktree and baseline before execution; these are references, not
instructions to reset or switch branches.

Scope amendment, 2026-09-08: the completed Codespaces feasibility pass is retained
in [the findings](runtime/CODESPACES_FEASIBILITY.md). The user subsequently
superseded the Docker deferral: both Docker and Podman must pass end to end on
the owned Kubuntu host reached through `ts ssh kubuntu`. Scoped installation and
host configuration there are authorized, with interactive sudo authentication
by the user when needed. Preserve unrelated remote workloads and local changes;
do not install Docker on the local Podman-only machine or weaken isolation.

## 1. Goal and finish line

Deliver a trusted-user LCM playground that restores the scientific showcase as
its public entrance, keeps the reusable toolkit in a developer gallery, and
launches registered presentations and workbenches against explicitly selected,
prepared, independently supervised Julia resources.

Finish with one migrated scientific showcase deck and one small scientific
workbench sharing real views, worker controls, diagnostics, and a private Julia
terminal. Demonstrate that losing a solver, terminal, UI host, or broker does
not take down the public site or an unrelated application run.

This is the completed execution contract; the ledger records acceptance evidence
and its limits. Pursuing it authorizes in-repository implementation and isolated tests;
it does not authorize public deployment, changes to unrelated services,
destruction of existing containers/data, or installing Docker on this machine.

### Locked release scope

- One coordinator and one local SQLite database; existing NATS/JetStream and
  artifact infrastructure remain in use. No PostgreSQL or Kubernetes.
- Linux first. Native Julia for trusted development/scientific execution;
  Podman and Docker adapters for approved container profiles.
- UI hosts run alongside the gateway in separate processes in v1. Scientific
  workers and terminal supervisors may run locally or on another computer.
  Remote workers connect outbound to the broker; no inbound Julia port or SSH
  orchestrator is introduced.
- Public static showcase, deck catalogue, workbench catalogue, and developer
  reference pages. Resource allocation and private runtime content require an
  authenticated principal.
- Operator-provisioned trusted users, not open signup or hostile public code
  execution. Containers are an operational isolation boundary, not a promise
  that arbitrary malicious code is safe on the host.
- Persistent worker registrations and run/lease bookkeeping. Full project
  CRUD, collaboration, autosave, billing, and persistent REPL memory are deferred.
  Application definitions and wire records leave room for a future project ID.
- One curated showcase deck adapted from the six ICHQP2026 sections, retaining
  their scientific narrative; no requirement to port every legacy control.
  Its live cases are lightweight visual interaction, line-parameter evaluation,
  and the existing OHL/UGC power-flow case. The terminal is also demonstrated.
- One `CableStudy` workbench uses those same view/operation adapters. The
  existing starter and hostile decks remain developer examples and additional
  catalogue registrations; no wholesale migration of every old presentation.
- Keep the current `lcm <feature> <action>` entry point. Add runtime actions
  there; do not undertake the separate top-level CLI reorganization.

### Explicitly outside this goal

New physics, a new layout engine, another theme system, automatic Makie PDF
snapshots, arbitrary remote machine provisioning, browser-supplied container
images/commands/mounts, remote UI-host placement, cross-host live process
migration, production HA, and public untrusted REPL access.

## 2. Reuse the existing foundation

Preserve and extend, rather than replace:

| Existing authority | Changes needed for v1 |
|---|---|
| `protocol/src/` | Versioned run, worker incarnation, assignment, preparation, and terminal-control records |
| `src/broker/` | Runtime client views, enforced targeting, reconnect/reconciliation |
| `worker/src/Consumer.jl` and `Executor.jl` | Profile-specific supervision, warm-state tracking, lease fencing |
| `worker/src/operations/` | Adapt existing scientific operations; separate their dependency loading |
| `src/workbench/` and `src/toolkit/` | Compose new portable controls; preserve existing shell contracts |
| `src/diagnostics/` | Register new owned components and test metadata/redaction/preview contracts |
| `_extensions/lcm-deck/` | Pass application-run context to live frames without changing Reveal layout ownership |
| `assets/brand.css` and `control-contract.css` | Remain the only shared palette and native-control state authorities |
| `src/container_runtime.jl`, `deploy/` | Reuse runtime detection, TLS profiles, artifact access, and lifecycle commands |

At baseline, the heartbeat advertised capabilities and slots, not executor warmth.
The baseline durable job request had no application run or enforced worker
assignment. Existing process supervision was useful but did not make a native
REPL a sandbox. These are actual gaps, not missing labels in the UI.

Recover narrative, assets, scientific inputs and useful preparation UX from
the old branch. Do not copy its global resource Observables, page machinery,
in-process scientific callbacks, or Serialization/Core.eval REPL transport.

## 3. Ownership and extension contracts

### Processes and deployment ownership

```text
Browser
  └─ Same-origin gateway: static publishing, authorization, live-route proxy
       ├─ Runtime coordinator: inventory, admission, leases, reconciliation
       ├─ UI host for application run A: Bonito sessions and view state
       └─ UI host for application run B: separate Bonito sessions and view state

Runtime coordinator / authorized runtime clients
  └─ NATS: bounded, versioned control and job messages
       ├─ Host agent → approved scientific executor profiles
       └─ Host agent → private, resource-limited Julia terminal containers
```

The gateway, coordinator and host agent never evaluate user code or load
scientific packages. UI hosts may load rendering packages, but not solvers.
A host agent is a restricted process supervisor, not a general remote shell.
Only it has permission to invoke its local container runtime; neither the web
publisher nor a user-code container receives a container socket.

Each dynamic application run receives one isolated UI host by default; all
live frames in its deck attach to that run. Static pages allocate none.
Admission limits prevent an arbitrary number of runs from exhausting the host.
Threads may accelerate an executor; they do not replace process isolation.

In v1 the gateway proxies local UI-host HTTP, WebSocket and asset paths using a
run-specific prefix. Remote scientific/terminal agents use outbound broker
connections. Do not tunnel Bonito/WebGL asset traffic through durable job
subjects.

### Small vocabulary, explicit state ownership

- **Application definition:** stable ID, title, kind, entry point, version,
  public/developer visibility, and named runtime requirements. Registration
  performs validation, not model loading or scientific execution.
- **Application run:** authenticated owner, application version, selected
  profiles, UI host, current assignments, lifecycle. A browser session and a
  persistent project are not interchangeable with a run.
- **Worker profile:** approved environment/image digest, operation capabilities,
  preparation recipe, resource budget and isolation policy.
- **Worker registration:** approved identity, allowed profiles and credential
  references. Discovery cannot approve a worker.
- **Worker presence:** current boot ID, liveness, capacity and per-executor
  preparation state. Persisted history cannot establish current readiness.
- **Lease:** owner/run/role, worker boot, assignment generation, bounded expiry,
  renewal and release. Executors reject stale or foreign assignments.
- **Job:** registered scientific operation, passive validated input and
  correlatable result. It is never an arbitrary Julia expression.
- **Terminal session:** private disposable process, separately authorized byte
  stream, explicit interrupt/restart and cleanup lifecycle.

Use immutable Julia records and concrete action types for these boundaries.
Keep current workbench `initialize`, `compose`, and `handle!` dispatch.
Add required profile hooks for validation, preparation, execution and cleanup
only where extension is actually needed; use RequiredInterfaces directly in
the lightweight layer where appropriate, without importing the engine for its
macros. Validate registrations at startup, including missing hooks, duplicate
IDs, incompatible versions and unsupported isolation requirements.

The owned root orchestrator enforces:

```text
authorize → validate → admit → acquire assignment → prepare
          → execute → persist terminal outcome → release/reconcile
```

Application hooks cannot replace authorization, lease checks or cleanup.
Use ordinary dispatch, not a plugin DSL, dynamic source evaluation, or a large
inheritance hierarchy. Wire strings map to registered types; never deserialize
arbitrary Julia types.

### Proposed source ownership

These are intended additions, not existing APIs:

- `src/applications/`: catalogue definitions, launchers and gateway/UI-host
  integration; reusable views remain outside page composition.
- `runtime/`: separate lightweight Julia environment for coordinator and host
  agent, SQLite persistence, authorization policy and supervision.
- `protocol/src/`: shared passive wire contracts only.
- `worker/`: reusable scientific worker core and profile-specific adapters/
  environments. Merely adding a profile flag while importing every engine
  would not meet the dependency boundary.
- `src/widgets/`: worker selector, preparation/status, diagnostics and terminal
  components, each with owned CSS and X-ray metadata.
- `presentations/showcase.qmd`, `src/workbenches/CableStudy.jl`: concrete
  consumers; no copied widget implementations.
- `dev/`: canonical developer landing/navigation; existing specimen URLs keep
  compatible aliases or links.
- `test/integration/`: one aggregate runtime-platform acceptance runner plus
  focused lifecycle, security, browser and runtime-adapter tests.

Exact file splits may change to fit the existing code. The dependency and
ownership boundaries may not be silently relaxed.

## 4. Operational rules to implement before scientific migration

### Identity and access

Provide a small, maintained reverse-proxy authentication deployment, using
Caddy HTTPS/basic-auth with provisioned accounts for the private v1 profile.
No custom password database or browser login framework is required.
[Caddy documents hashed credentials and the authenticated user identity](https://caddyserver.com/docs/caddyfile/directives/basic_auth).

Only a trusted, non-public proxy connection may assert that identity; strip
client-supplied identity headers, and make direct gateway bypass impossible in
the deployed profile. Separate public static paths from protected application,
runtime API, WebSocket, terminal and private artifact paths. Authorize each
operation against principal, run and lease, not just the initial HTML request.
Audit existing Bonito job routes and v1 submission paths as well: leaving an
old unauthenticated calculation callback reachable would bypass the new policy.

Add explicit origin and CSRF protection to mutations and WebSocket upgrades.
Use fixed configured public origins, never a reflected client Host header.
Credentials, terminal input, secrets and raw private paths are excluded from
X-ray and ordinary logs. A development identity adapter is explicitly
loopback-only and opt-in; it is not a deployment authentication mechanism.
A future institutional/OIDC adapter may replace the proxy identity source
without changing ownership checks.

### Inventory, placement and leases

SQLite is authoritative for approved registrations and allocation bookkeeping;
NATS is transport and existing durable job/result storage. Do not introduce a
second independently authoritative lease table in JetStream.

Registration binds the broker-authenticated worker identity to permitted
profiles and subjects. Payload `worker_id` alone is not authentication.
Use least-privilege per-agent credentials or an equivalently enforced mapping;
agents cannot announce as another worker or publish grants to themselves.

Selectors support automatic, pinned-worker and dedicated-run placement.
Capacity reservations are transactional. A lease is usable only after the
target acknowledges it. Control actions have request IDs and idempotent
responses; grant expiry/renewal uses bounded durations and local monotonic
timers rather than trusting browser clocks.

Use explicit dimensions instead of one misleading green status:

- liveness: unknown / online / stale / offline;
- preparation: cold / preparing / ready / failed, per profile and executor;
- allocation: free / reserved / busy / draining.

No matching worker is a finite unavailable state with a reason. Pinned-worker
failure never silently chooses another worker. Revocation/draining rejects new
work; terminate-existing-work is a separate explicit action.

Lease expiry prevents new starts, cancels or stops owned execution within its
configured bound, and fences late results. Renewals during broker loss are
not assumed to succeed. A coordinator restart reconciles fresh agent reports;
a saved lease row is not proof of a live process.

### Versioned job delivery

Introduce a v2 contract for assigned jobs/control, with separate subjects and
consumer filters. Preserve v1 consumers during migration; mark legacy workers
as unable to provide enforced selection until upgraded.

A granted assignment binds run, role, worker boot, executor generation and
operation/profile compatibility. Only its intended worker/assignment receives
the targeted job. Do not pull from a shared queue and discard another worker's
job, or rely on a selector's label for placement.

Keep durable result-before-ack and bounded retries. Scientific operation
retries must be idempotent or explicitly disallowed after uncertain execution;
do not promise exactly-once external effects. Cancellation, cache keys, result
authorization and stale-result rejection carry the same ownership context.
Never blindly replay uncertain terminal input.

### Preparation, caches and cleanup

Separate package artifacts, warm processes, prepared models, numerical results
and browser view state. Their lifetimes and cache keys are different.

Readiness requires successful preparation of the actual assigned executor,
with environment fingerprint, profile version, representative workload and
input identity recorded. Worker replacement invalidates its warm status.
A representative warmup reduces expected compilation; it is not a guarantee
that every future Julia specialization is compiled.

Deduplicate concurrent preparation for the same compatible ownership/key.
Bound queue lengths, preparation time, cache bytes and cache lifetime. Do not
share mutable prepared models across unrelated runs. Display stages, elapsed
time, failure and retry without blocking Bonito callbacks.

No scientific preparation on page construction or slide entry. The user
explicitly prepares/launches the run. Last successful data remains visible
with provenance while another request is pending.

All resource-owning paths have idempotent cleanup: normal stop, rejected
launch, failed preparation, cancellation, disconnect expiry, agent restart and
coordinator recovery. Agents track only their owned process/container IDs and
scratch directories. Reconciliation never removes unrelated containers.

### Terminal boundary

Use pinned, locally served xterm.js plus an actual Julia REPL attached to a PTY;
do not build another Julia parser/line editor. Browser keys and terminal size
travel on a separately authorized, bounded session channel through the gateway
and agent, not the durable scientific job stream. Remote relay subjects have
per-session permissions, ordering/sequence checks, bounded chunks and explicit
gap/disconnect behavior. Control traffic must not wait behind terminal output.

Selecting a worker/profile creates a private terminal container there; it
does not attach to a shared scientific executor. Default policy: no network,
no broker credentials, no home-directory mount, read-only environment,
quota-limited scratch, unprivileged user, dropped capabilities and enforced
CPU, memory and PID limits. Prebuilt profiles provide required packages.

Preflight must verify that the host actually enforces the required limits.
Rootless Podman controller availability depends on host configuration; if a
required limit is unavailable, shared terminal launch fails closed, while the
rest of the application remains usable.
[Podman documents the available isolation and resource controls](https://docs.podman.io/en/latest/markdown/podman-run.1.html).

Allow one writer per terminal; reconnect only as the same authorized owner.
The terminal preserves state only while its process survives. Restart creates
a clean process; recovery never pretends that lost variables were restored.
Use bounded interrupt followed by explicit hard-stop/restart escalation.
Disconnect grace, idle TTL, output rate/scrollback and scratch quota are
configured and tested. Do not record raw keystrokes in durable job/audit logs.

Render terminal output as untrusted text, with unsafe clipboard/control
extensions disabled and links restricted. Key handling must stop presentation
shortcuts while the terminal owns focus. Terminal code and secrets stay outside
X-ray.
[xterm.js explicitly requires secured terminal transport and careful handling of output](https://xtermjs.org/docs/guides/security/).

## 5. Ordered implementation work packages

Each package includes tests and a usable checkpoint. Do not defer all hardening
to the last package. Record its evidence before marking its checklist complete.

### P0 — Baseline and acceptance inventory

- [x] Record current branch/worktree changes, manifests, tool versions, test
      commands, available Podman capabilities, and Docker-runner availability.
- [x] Run existing protocol, publisher, worker, presentation, ribbon and X-ray
      gates. Record pre-existing failures separately from new regressions.
- [x] Inventory current registered components and legacy showcase cases. Freeze
      the migration checklist: six narrative sections, three live case families,
      one workbench and one terminal implementation.
- [x] Create the runtime requirement/test matrix and versioned config schema.

**Gate:** reproducible baseline, exact case inventory, and no unexplained
changes to user-owned files or running services.

### P1 — Contracts, identity and an isolated UI-host vertical slice

- [x] Add application/profile registration contracts, principal/run identity,
      fixed orchestration hooks and SQLite migrations/transaction tests.
- [x] Add trusted-proxy and explicit local-development identity adapters;
      prove authorization with two distinct test principals.
- [x] Launch a mock run's UI host through an allowlisted local supervisor.
- [x] Proxy its HTTP, WebSocket, assets, theme messages and X-ray through the
      same-origin run namespace. Route two live deck frames to the same run.
- [x] Bound startup and concurrent host count; implement stop and orphan
      reconciliation. Preserve a lightweight unavailable/restart surface.

**Gate:** two application runs cannot reach each other's sessions. Killing one
UI host leaves the gateway and other run usable. Reconnect to a surviving host
preserves its session; host death is explicitly reported as lost volatile state.
No scientific packages or broker availability are needed for this test.

This gate resolves the riskiest integration before porting pages. If Bonito
cannot satisfy the namespaced proxy/asset contract without violating owned
boundaries, document the evidence and revise that adapter before proceeding.

### P2 — Public showcase, developer gallery and application catalogue

- [x] Restore a public scientific landing page using the old showcase's intent
      and existing branding, not its runtime implementation.
- [x] Make `/presentations/` and `/workbenches/` real registered-application
      launchers. Definitions declare requirements without loading their engines.
- [x] Move the foundation/frontface to `/dev/`; keep developer access public
      but separate from the primary navigation.
- [x] Preserve existing template/widget/deck URLs and Source links through
      explicit compatible routes; validate generated links after Quarto build.
- [x] Register the showcase skeleton, template workbench and existing developer
      decks. Keep selection correct when additional decks are registered.
- [x] Show unavailable runtime capabilities honestly, with no implicit startup.

**Gate:** the complete public/developer navigation and static slides work with
NATS and every worker stopped. Both themes, document scrolling, print
placeholders, fragment lists, pointer behavior and live-frame identity retain
their existing contracts.

### P3 — Worker control plane, targeting and diagnostic components

- [x] Implement approved registration/enrollment, authenticated announcement,
      boot identities, TTL pruning and inventory reconciliation.
- [x] Implement v2 lease grants/acknowledgement/renewal/expiry, admission,
      dedicated and pinned placement, and targeted job delivery.
- [x] Add control panel actions to approve/disable registrations, inspect
      capacity, assign roles, drain and release; credentials are references,
      never browser-visible broker secrets.
- [x] Provide launch configuration for approved native/Podman/Docker agents.
      Host installation and image provisioning remain operator actions.
- [x] Build reusable WorkerSelector, PreparationStatus and WorkerDiagnostics
      components against the same runtime client, not per-page discovery code.
- [x] Correlate events by run, worker boot, lease, executor, job and stage.
      Bound/redact log buffers; mark dropped events and reconnect gaps.
- [x] Keep terminal/log floods independent of heartbeat and control scheduling.

**Gate:** duplicate enrollment, spoofed identity, incompatible profiles,
oversubscription, stale grants, cross-owner actions and wrong-worker execution
are rejected. Pinned-worker loss is visible without silent reassignment.
Duplicate/delayed messages and coordinator restart do not double-allocate.
Inventory/control behavior is identical from the panel, gallery and workbench.

### P4 — Profile supervision and honest preparation

- [x] Factor scientific dependencies into explicit approved profiles, reusing
      current operation adapters. Separate line parameters from power flow.
- [x] Keep UI rendering and lightweight local animation separate from those
      processes; move mathematical animation computation to an executor if it
      can block the UI host.
- [x] Extend the existing executor supervisor for assignment generations,
      readiness/progress, resource policy and bounded cancellation/cleanup.
- [x] Add allowlisted terminal-container launch/PTY lifecycle to the host agent,
      with limit preflight; no browser UI is required for its first tests.
- [x] Implement preparation deduplication, cache ownership/invalidation,
      deadlines and finite error/retry states.
- [x] Add per-profile resource budgets and per-run/user admission limits.

**Gate:** a genuinely prepared process reports ready; its replacement reports
cold. Power-flow failure/preparation cannot block line-parameter execution,
heartbeats or public navigation. Repeated failed starts leave no owned stale
processes, containers, leases or scratch files. Unsupported sandbox limits
disable only that capability, not the site.

### P5 — Portable live Julia terminal

- [x] Compose the terminal component with worker selection, readiness,
      connect/disconnect, interrupt, restart, clear and status.
- [x] Implement the authorized local and remote stream relay with backpressure,
      bounded reconnect and no durable input replay.
- [x] Demonstrate the same component in the developer gallery, a deck and a
      workbench; no gallery-only styles or page-specific terminal behavior.
- [x] Exercise multiline Julia input, history within the session, completion,
      Unicode, resize, terminal exit, interrupt and forced replacement.
- [x] Test output floods, allocation/PID limits, unauthorized attachment,
      expired ownership, network/mount restrictions and teardown.
- [x] Verify light/dark switching, compact geometry, keyboard ownership and
      X-ray redaction.

**Gate:** two owners receive separate terminal state. Infinite work, process
exit, output flood or container death cannot freeze the UI or another worker.
All three hosts use the identical component, and presentation navigation does
not consume terminal keystrokes.

This is the first complete operational demonstration, before legacy scientific
UI migration.

### P6 — One real deck and one real workbench

- [x] Adapt the six ICHQP2026 narrative sections into Quarto `showcase.qmd`;
      retain approved assets and scientific meaning. Register it as public.
- [x] Extract portable live views for lightweight interaction, line parameters
      and the OHL/UGC power-flow case. Reuse existing validated numerical
      operations and artifact retrieval rather than porting old global caches.
- [x] Declare separate runtime roles/profiles, expose selection and explicit
      preparation, and show per-role readiness/progress/failure.
- [x] Compose `CableStudy` from the same views, forms, job controls, terminal
      and diagnostics. No duplicated scientific operation or control markup.
- [x] Keep result provenance and last-good data visible; discard stale
      completion after changed inputs, lease generation or run replacement.
- [ ] Provide explicit versioned input export/import for the small workbench
      if needed for the demo; this is not a new project-management service.

The curated v1 specimen does not currently need file import/export: its small
input set is explicit and repeatable in both consumers. This optional item does
not imply project persistence. The P6 numerical/failure gate below has passed
with finite native test executors. P7 subsequently passed both physical engine
matrices and the separate-computer consumer gates. Results and limitations are
recorded in the ledger.

**Gate:** perform the same real calculation from both hosts and compare its
normalized inputs/results with the existing worker tests within their declared
tolerances. Repeat cold/warm preparation, cancellation/retry and worker-loss
cases. Preparation occurs before the live demo, not as an unexplained pause on
slide entry. Static deck navigation and PDF placeholders remain independent.

### P7 — Release verification and handoff

- [x] Add one aggregate `run-runtime-platform.sh` acceptance entry point with
      isolated namespaces, temporary storage and cleanup traps.
- [x] Derive the component conformance inventory from registrations; require
      every new owned component to supply its normal render path, metadata,
      CSS ownership and test fixture. No handwritten gallery replicas.
- [x] Run registered components in gallery/workbench/deck contexts where
      supported, both themes, hover/focus/selected/disabled states and repeated
      mount/unmount. Include X-ray preview/reset and metadata redaction.
- [x] Run all existing protocol, worker, artifact, shutdown, presentation,
      ribbon and X-ray regression gates. Update ARCHITECTURE and VERIFICATION.
- [x] Run native trusted-development checks, actual Podman checks, and actual
      Docker checks on an available Docker host/runner. Do not install Docker
      on the Podman-only machine or call adapter mocks an integration pass.
- [x] Exercise remote topology with separate process/container identities over
      authenticated TLS, outbound-only agents and no shared filesystem.
      A physical second-computer rehearsal is a separately recorded deployment
      check; the portable transport contract must already be automated.
- [x] Deliver exact CLI/config examples, environment/profile inventory,
      backup/restore and migration instructions, operator troubleshooting,
      measured cold/warm times and resource requirements.
- [x] Record every acceptance result, test artifact, known limitation and any
      unavailable external verification.

**Gate:** the final scenario below passes and all mandatory release gates have
evidence. Missing Docker access is reported as pending cross-runtime
verification, never silently converted to a pass or a fully verified release.

## 6. Final acceptance scenario

Run this as an automated scenario plus a short manual test drive:

1. With NATS stopped, open the public landing, browse developer references,
   choose a registered deck and navigate/print its static content.
2. Start the approved runtime stack. Authenticate two principals; enroll a local
   agent and a remote-topology agent. No unauthorized self-enrollment succeeds.
3. Open a presentation run and a workbench run. Select their profiles. Prepare
   the required executors, observing actual stages and independent readiness.
4. Run line parameters and the power-flow case; verify numerical regressions
   and result provenance. Reuse a compatible warm preparation and measure it.
5. Open a private Julia terminal in each host. Verify completion, input,
   resizing, independent state and focus-safe slide navigation.
6. Kill a power-flow executor during work. Its status fails explicitly; the
   line-parameter role, terminal, menus and public site remain responsive.
7. Flood/terminate one terminal within its enforced bounds. Recover with a clean
   restart; the other terminal and all scientific workers remain usable.
8. Kill one UI host. Its run reports loss/restart; the other run and public
   publisher remain usable. Do not claim recovery of unsaved session memory.
9. Disconnect NATS long enough to expire leases, then restore it and restart
   the coordinator. Reconcile without duplicate jobs, stale authority, falsely
   warm replacement processes or cross-owner results.
10. Close/release the runs and wait the configured disconnect grace. Compare
    owned-resource inventory against baseline: no new orphan process/container,
    lease, temporary upload or scratch directory remains.
11. Repeat visual interactions in both themes, including X-ray reset,
    presentation overview, fragment lists and linked-placeholder PDF export.

Tests use declared bounds, not indefinite waiting. The baseline should retain
the existing 2-second heartbeat / 10-second presence timeout unless measurements
justify a documented change. Set profile-specific preparation/job deadlines and
measure browser responsiveness under load; numerical execution duration is not
allowed to determine whether menus or status controls respond.

## 7. Execution discipline, compatibility and completion reporting

- Implement in P0–P7 order, allowing only independent work within a passed
  boundary. Each work package lands with its tests; do not mark scaffolding as
  a completed capability.
- Extend current architecture tests to forbid scientific dependencies in the
  gateway/coordinator, direct NATS credentials in browser output, arbitrary
  evaluation on scientific subjects and application-local palette copies.
- Use versioned, backed-up SQLite migrations. Add v2 NATS subjects/streams
  without purging v1 jobs/results; document draining and rollback explicitly.
  Refuse unsupported protocol/schema versions instead of guessing.
- Feature-gate new launchers until their gates pass; static aliases remain
  available. Never restore the legacy runtime as an implicit failure fallback.
- Do not add new domain features while pursuing the platform goal. Unexpected
  physics changes, public hosting, stronger sandbox infrastructure, incompatible
  proxy requirements or unavailable mandatory host capabilities require a
  scoped decision and evidence, not an unannounced architectural compromise.
- Maintain an implementation ledger with checkbox state, commands, results and
  blockers. This plan is not a substitute for a verification report.
- Final handoff distinguishes implemented, verified locally, verified on another
  runtime, and deferred. No claim of production/public sandbox readiness.
- Do not automatically commit/push, deploy publicly, or remove pre-existing
  resources. Use isolated labelled test resources and clean up only those owned
  by this execution.

### Goal text for a subsequent execution request

> Implement Runtime Platform v1 according to
> playground/RUNTIME_PLATFORM_PLAN.md, completing P0–P7 with the specified
> ownership, security, visual and recovery contracts. Deliver the restored
> scientific showcase, registered deck/workbench launchers, runtime control
> panel, enforced worker selection and preparation, reusable diagnostics and
> private Julia terminal, one real deck and one real workbench, plus repeatable
> tests and operator documentation. Preserve existing behavior and user changes.
> Report unavailable external verification explicitly; do not broaden scope or
> weaken mandatory isolation to obtain a nominal pass.
