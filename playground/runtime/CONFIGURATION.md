# Runtime configuration v1

The tables below, remote scientific preparation and private terminal lifecycle
are implemented. Configured private terminal transport is advertised separately
from executor eligibility: only an approved installed image on a verified managed
host can be used. Consult `RUNTIME_PLATFORM_PROGRESS.md` for measured acceptance;
development transport tests alone do not certify physical isolation or deployment.

Use TOML with integer `schema_version = 1`. Reject unknown keys, duplicate
identities, invalid units/ranges and unsupported versions. Paths resolve
against the configuration file directory, never a browser-provided directory.
Keep secrets in private files or environment references, outside source control.

| Table | Required fields / defaults | Validation |
|---|---|---|
| root | `schema_version = 1`, `enabled = false` | Unsupported versions/unknown keys fail |
| gateway | `listen_host = "127.0.0.1"`, `port = 8080`, `public_origin` | Exact configured origin; no path/query/credentials; HTTPS except explicit loopback development |
| identity | `mode = "proxy"`, `proxy_peers`, `proxy_key_file`, `administrators` | Trusted peer AND private proxy key AND asserted principal; no client-controlled identity |
| storage | `database`, `scratch_root` | Private owned directories; SQLite schema version verified |
| publisher | optional `site_directory` | Existing Quarto build; no symlinks or reserved runtime routes |
| limits | `max_runs = 8`, `max_runs_per_owner = 2`, `startup_seconds = 60`, `shutdown_seconds = 5`, `disconnect_grace_seconds = 60` | Positive finite bounded values; admission before spawn |
| control | optional `config_file` | Separate strict server-owned profile/trust/broker configuration; no connection while checking |

Development mode is `identity.mode = "local-development"`, explicitly chosen,
with a configured local principal and a loopback-only listener/peer/origin.
It does not accept identity headers and cannot be enabled on a public listener.

Profiles are separately registered, immutable approved records: profile ID,
version, environment fingerprint/image digest, capabilities, launch kind,
resource limits, prepare/job deadlines and cleanup policy. Browsers select
only registered IDs; never commands, mounts, paths or image URLs.

Application requirements name roles and compatible profiles, not machine
addresses. Placement selects automatic, pinned or dedicated assignment.
Run ownership and a future project reference are distinct.

Unknown tables are rejected, not silently ignored. In particular, terminal
launch is not enabled by inventing a `terminal` table.

P2 implements `publisher` and the `lcm runtime check/start/status` actions.
For an explicit loopback-only test drive, from `playground/`:

```sh
./lcm playground build
./lcm runtime check --config runtime/local.example.toml
./lcm runtime start --config runtime/local.example.toml --xray
```

This opt-in path serves the static build and run-prefixed owned UIs without
importing Bonito into the gateway. The old `lcm playground start` remains the
developer publisher; it is not the authenticated runtime deployment entry.
The local example is not a public-deployment identity policy. Its state is
written under the ignored `runtime/state/` directory.

Build/install the UI environment before admitting application runs. Its package
image now includes representative disconnected UI rendering, so common form and
workbench compilation does not fall on the first browser request. Package build
can take longer after a source/dependency change; it does not contact a broker,
launch workers or prepare numerical models. Actual host startup still has its
configured finite readiness deadline.

For the approved trusted-proxy profile, use the paired templates and explicit
validation/rehearsal checklist in [PROXY.md](PROXY.md). They do not install a
proxy, expose a listener publicly or certify a deployment.

## Worker control

Opt in from the runtime TOML with:

```toml
[control]
config_file = "control.toml"
```

That separate file is also schema version 1. It contains:

- `[broker]`: `url`, `password_file`; optional `ca_file`,
  `certificate_file` + `key_file`, `server_name`, and explicit
  `allow_loopback_plaintext = true` only for literal-loopback development.
  Remote endpoints require TLS. Secrets must be nonempty regular files, not
  symlinks, with mode 0600 or stricter. URLs cannot contain credentials.
- `[[profiles]]`: `id`, `environment`, lowercase SHA-256 `fingerprint`;
  defaults `version = "1.0.0"`, `kind = "scientific"`,
  `isolation = "trusted_process"`, `preparation = "default"`,
  `protocol_version = 2`. Scientific profiles require distinct `operations`.
  A container environment must be digest-pinned and match its fingerprint.
  Terminal profiles require container isolation and no scientific operations.
  Optional `[profiles.budget]`: `cpus`, `memory_bytes`, `pids`,
  `scratch_bytes`, `prepare_seconds`, `job_seconds`.
- `[[workers]]`: `id`, opaque `credential_ref`, allowed profile ID array
  `profiles`, and `capacity = 1` by default. Up to 128 distinct provisioned
  workers and 64 profiles. Each credential reference belongs to one worker.
- `[assignments]`: optional `total = 128`, `per_owner = 8`, `per_run = 4`.
  Pending grants and cleanup/reconciliation still occupy capacity.

All paths are relative to the control file. Native environment references are
passive: the coordinator neither imports them nor claims to have verified their
installed contents. Actual profile preflight/preparation is an agent duty.

From `playground/`, after creating the operator-owned files:

```sh
./lcm runtime check --config runtime/local.example.toml
./lcm runtime permissions --config runtime/local.example.toml
# Install the printed users and distinct matching passwords in the broker.
# Then explicitly create/verify only the bounded v2 worker streams:
./lcm runtime provision --config runtime/local.example.toml
./lcm runtime start --config runtime/local.example.toml
```

`permissions` is read-only and prints environment-variable references, never
secret contents. `provision` connects and creates missing per-worker v2
streams/consumers; existing definitions must match. It never purges v1 streams,
silently resizes existing streams, or modifies broker users. A partial failure
can leave successfully created v2 streams; rerun after correcting configuration.
Host installation, broker credentials and image provisioning remain operator
actions, not browser features.

Configuration only makes an identity eligible for enrollment. The operator
must separately enroll it (pending) and approve its exact trust binding. A
changed persisted trust binding is rejected at startup rather than silently
expanding permissions. Removing an identity from the configuration stops new
discovery; historical registrations/leases remain for reconciliation.

The running coordinator uses a 2-second probe interval, 10-second presence
expiry, 10-second renewable lease and 2-second acknowledgement bound. Initial
connections are finite background attempts; an unavailable broker does not
prevent public HTTP startup. Approval, liveness, reservation and preparation
are separate dimensions. Advertisement never means ready.

Explicit service startup compiles the installed control/scientific orchestration
entry points before announcing availability. This performs no driver hooks,
assignments, model preparation or numerical work. Protocol codecs also declare
their passive precompilation paths. A cold Julia process can therefore spend time
compiling during startup; it must not charge that known first-use work to a live
lease acknowledgement. The 2-second acknowledgement bound is unchanged. Model
preparation remains a separate, explicit user action in the selected executor.
The configured HTTP entry is likewise compiled before its listener opens; the
CLI does this before starting broker scheduling. This is an inert compiler
request, not a synthetic authenticated request or a longer browser timeout.

### Protected API

All endpoints below use the same principal, origin and CSRF policy as owned
application runs. Mutations require `Origin` and `X-LCM-Request: 1`. Responses
exclude environment paths, broker credentials and another owner's run context.

| Endpoint | Action |
|---|---|
| `GET /runtime/api/control` | Shared profile/inventory dimensions; provisioned enrollment choices for administrators only |
| `GET /runtime/api/control/events?after=N&epoch=UUID` | Bounded owner-filtered structured events, incarnation/cursor and eviction/reconnect gap |
| `GET /runtime/api/workers` | Approved/pending registrations, liveness, occupied slots; no implied preparation |
| `POST /runtime/api/workers` | Administrator enrollment: `worker_id`, `request_id`; trust comes from configuration |
| `PATCH /runtime/api/workers/ID` | Administrator `state` (approved/draining/disabled), `expected_revision`, `request_id` |
| `GET /runtime/api/runs/UUID/assignments` | Owned assignment history |
| `POST /runtime/api/runs/UUID/assignments` | Owned `role`, `profile`, `placement`, `request_id`; reserve and request acknowledgement |
| `GET /runtime/api/assignments/UUID` | Owned current state; usable lease is not a prepared executor |
| `DELETE /runtime/api/assignments/UUID` | Owned explicit release with `request_id`; wait for cleanup acknowledgement |
| `GET /runtime/api/assignments/UUID/science` | Request/read bounded executor status; never start or prepare a child |
| `POST /runtime/api/assignments/UUID/science` | Explicit `prepare` or exact-request `cancel`, as described below |

Durable submission, result download and cancellation endpoints are specified in
[Protected scientific job API](#protected-scientific-job-api).

Placement is `{"mode":"automatic"}`,
`{"mode":"pinned","worker_id":"worker-a"}`, or
`{"mode":"dedicated","worker_id":null}` (optionally pin its ID).
Duplicate/unknown JSON fields and excessive nesting are rejected.
Disabling or draining a registration rejects new work but does not implicitly
kill existing execution; release is a separate operation. Broker credential
revocation is also a separate operator action.

The capability endpoint and private control inventory use one service declaration.
`preparation_control` and `assigned_execution` indicate configured, open service
owners, not worker readiness or isolation evidence. Run additionally requires an
acknowledged live assignment, an allowed operation and fresh preparation. Closing
the job service removes execution availability from both responses.
`private_terminal` likewise describes the configured relay, not an eligible or
prepared worker. The selected terminal image and actual managed host must pass
their separate admission and kernel-isolation checks.

### Explicit scientific preparation

Preparation uses the assigned profile's approved recipe and passive JSON inputs:

```json
{"action":"prepare","parameters":{},"request_id":"YOUR_REQUEST_UUID"}
```

The UUID must be generated once and retained for an explicit retry. Inputs are
bounded to 64 KiB. The gateway checks the owner, live acknowledged lease and
scientific profile before publication; the agent checks the same exact worker
boot, coordinator, assignment and generation before admitting work. A 202 HTTP
response means the request was recorded, not that the executor has prepared.
Conflicting work is rejected, not placed in an unbounded queue. Matching pending
preparations share the existing task.

Poll the GET endpoint for `channel`, `phase`, `preparation`, accepted/reason,
executor identity/generation, progress, elapsed seconds, output-line count and
fixed failure codes. Raw child output, exceptions and paths are not included.
Only an explicit retained-state inspection of the current child can report
`ready`. The coordinator subtracts the full broker round trip from the report's
finite validity (at most five seconds), and the shared browser client subtracts
HTTP latency again. Reading a report never renews it or warms the model.
An unexpired report may remain visible during a read-only successor query;
preparation/cancellation, lease loss, expiry or a transport failure invalidate
readiness. A failed status query never silently prepares a replacement child.

Cancellation names the exact `current_request_id` shown by the status response:

```json
{"action":"cancel","target_id":"CURRENT_REQUEST_UUID","request_id":"NEW_CANCEL_UUID"}
```

A delayed cancel cannot stop a successor. Cancellation keeps the assignment;
release is the separate lease operation. Mutations are never replayed
automatically after an uncertain reply. Each lease retains at most 256 mutation
ID/digest tombstones, not an unbounded history of preparation inputs; after that
limit, release and reacquire the role before issuing a new mutation. Status reads
do not consume that history. Exact retries remain inert within the live lease.

Preparation traffic has its own bounded NATS connection and scheduler, separate
from heartbeat/lease control and durable job transport. Update the broker users
with the output of `lcm runtime permissions` when upgrading: the coordinator
publishes `lcm.science.v2.*.command` and subscribes to reports; each worker receives
only its own command subject and publishes only its own report subject. An agent
cannot grant authority to itself. Changing the application does not install
these broker permissions automatically.

### Control panel and reusable views

Open `/runtime/control` after authenticating. It lists only runs visible to the
principal. Choose a run to see selectors for its registered roles, or follow
`/runtime/control?run=RUN_UUID` from its recovery page. Opening this page allocates
no worker. The administrator may enroll a provisioned identity as pending, then
approve, drain or disable its registration. No credentials or launch commands
are accepted from the page. The owned-run list is a page snapshot; reload it to
refresh the list. Worker inventory and diagnostics poll live.

Bonito applications use the exported components from the publisher package:

```julia
using LineCableModelsPlayground, UUIDs

client = RuntimeClient(UUID(ENV["LCM_RUN_ID"])) # supplied by the owned UI host
selector = WorkerSelector(client, :parameters; profiles=("line-parameters",))
status = PreparationStatus(client, :parameters; parameters=Dict())
diagnostics = WorkerDiagnostics(client)
panel = WorkerControlPanel(client, selector)
```

Compose the individual values or the complete panel in an ordinary Bonito DOM
or workbench view. Presentation live frames consume the same components; no
deck-specific selector exists. `RuntimeClient()` is inventory-only and disables
assignment. The application must independently declare its roles and approved
profiles in its registered definition; changing a widget's choices cannot grant
execution permission. The gallery uses inventory-only selectors and allocates no
scientific resources.

All surfaces use `assets/runtime-client.js` and `assets/runtime-controls.js`, with
layout in `assets/runtime-controls.css` and palette/native-control states from
the existing brand and form styles. Bonito wrappers add source-scoped X-ray
metadata, not credentials or event payloads. Browser clients coalesce reads per
run/window, retain drafts and last-known stale status, separate bounded event
and preparation polling from inventory, and detach when their final component is
removed. Hidden pages and stale inventory clear browser readiness immediately;
finite readiness timers also expire without waiting for another poll.

An uncertain action is never retried automatically. **Retry same action** reuses
its original request ID and input; **Keep current state** dismisses that pending
UI choice without cancelling any action already accepted by the coordinator.
Inspect the refreshed assignment before changing it. Offline pinned workers stay
selected and visibly unavailable; there is no silent placement fallback.

Preparation displays **unknown** until fresh executor evidence arrives, then
**cold**, **preparing**, **ready**, or **failed** as appropriate. The component
offers explicit **Prepare executor** and **Cancel preparation** actions using
the same client and renderer in the panel, Bonito views and presentation frames.
Its configured input values are excluded from X-ray; command names and types
remain inspectable. A heartbeat, installed profile, or acknowledged lease never
produces a false ready indicator.
Structured diagnostics contain registered control-event codes and authorized
identities only, not arbitrary executor output or terminal history.

### Explicit calculations and retained views

`ScientificJob` composes with those same controls. Create it inside the owning
Bonito session, render the control once, and project its result into persistent
plots/tables without replacing their surrounding DOM:

```julia
parameters = Observable(Dict{String,Any}("value" => 3))
calculation = ScientificJob(client, :parameters, "system.echo"; parameters)
display_text = map(session, calculation.result) do result
    result === nothing && return "No successful result"
    suffix = result.current ? "current" : "previous result · outdated"
    return string(result.value, " · ", suffix)
end
DOM.div(calculation, DOM.output(display_text))
```

The example illustrates composition; it does not register the operation or
grant a role access to it. Use the operation and passive inputs approved by the
concrete application's definition. The configured `assigned_execution` capability
describes service availability, not worker preparation or deployment certification.

When fields drive the input Observable, mark their shared owner and field region:

```julia
DOM.div(calculation,
    DOM.div(input_widgets...; var"data-runtime-input-fields"=""),
    DOM.output(display_text);
    var"data-runtime-input-scope"="")
```

`input_widgets` are the instance's ordinary bound controls. The shared job binding
observes their native input/change events, disables Run immediately and waits for
the latest ordered Bonito acknowledgement containing the canonical input draft.
Older and duplicate acknowledgements cannot unlock or replace a newer draft.
Keep display-only selectors outside the marked field region. Each scope owns one
job; nested scopes retain their own input ownership. No arbitrary wait duration,
duplicated JavaScript domain validation, global event interception or implicit Run
is involved.

- Mounting, changing input values, preparing a worker, revisiting a slide, and
  switching theme never submit a job. **Run calculation** captures explicit intent.
- `result` holds `nothing` or `ScientificResult(provenance, value, current)`.
  The existing assigned-result record supplies job, input hash, worker/boot,
  assignment generation, prepared executor/model and schema identity. The control
  shows that provenance and a bounded JSON preview; the Observable retains the
  complete passive result. Inputs and values are excluded from X-ray metadata.
- Input changes retain the last successful value but mark it outdated. Readiness
  expiry, offline status or lost evidence also prevents a claim of currency.
  Known input/assignment/executor replacement and explicit revocation prevent
  that job's late completion from replacing the view. Identical input echoes
  are not new drafts. New successful explicit work can replace historical data.
- **Retry acknowledgement** reuses the captured request ID and original inputs,
  including when fields have since changed. An unconfirmed submission blocks a
  new Run. **Cancel job** shows requested and acknowledged cancellation separately
  from the durable terminal outcome. **Refresh job** reads without resubmitting.
- Successful results use the existing inline/private-artifact path, never a URL
  supplied by the result. Failed reads retain the prior value. The browser-to-Julia
  projection rechecks the passive wire shape, run/role/operation, current input
  hash and render-local draft token. It is display data, not an authorization
  token; the gateway remains the authority for execution and downloads.
- Closing/removing a control or departing frame stops its polling/subscriptions,
  not its computation and not sibling controls. Owned-run teardown and explicit
  cancellation remain separate. No result is automatically restored after reload.

### Agent configuration and ownership

`lcm runtime check-agent --config agent.toml` validates the separate agent
file without starting anything. Its schema-1 root accepts only `agent`,
`broker`, `profiles`, and optional `artifacts`. Broker and profile tables use the same grammar
above, but reference that machine's worker credential and local environments.
`[agent]` requires `worker_id`; optional `capacity = 1`,
`scratch_root = "state/agent"`, and
`container_runtime = "auto"` (`"podman"` or `"docker"` to pin selection).
Runtime selection here is only configuration, not a container preflight pass.

For a trusted native scientific profile, compute its approved source digest from
the target checkout (from `playground/`):

```sh
./lcm runtime fingerprint --project worker/profiles/line-parameters
./lcm runtime fingerprint --project worker/profiles/power-flow
```

This read-only action requires no broker configuration. It prints only the digest
and covers the selected manifest, mutable local package sources and declared
cross-package includes. It does not start Julia children or import either engine.
`check-agent` remains a passive configuration check; it does not imply that this
source verification, resource preflight or representative preparation succeeded.
See [scientific profile source and cleanup contracts](../worker/profiles/README.md)
for exact inputs, supported filenames, finite limits and the trusted-depot boundary.

`AgentService` owns announcement, lease expiry and acknowledgement scheduling.
It accepts one approved `AbstractAgentResources` implementation, whose required
hooks are `installed_profiles`, `bind_agent!`, `recover_owned!`, `release_owned!`, and
`close`. Hooks are checked before startup. Recovery must finish before the
first announcement. Cleanup tasks do not block heartbeat processing, and
unresolved cleanup keeps its slot occupied. A replacement coordinator is not
announced to until old owned resources have been cleaned. Repeated concurrent
shutdown calls join the full teardown; failed cleanup stays an explicit failure.
A subsequent close can finish the same unresolved cleanup after the driver recovers,
but can never reopen admission or reuse that closed agent incarnation.

An empty installed-profile report is allowed: the supervisor remains visible,
but cannot receive a compatible assignment. Installed profiles are not prepared
executors. The production schedulers have been verified over real TLS with
explicit no-executor resource fixtures. The shared managed scientific driver and
agent CLI are implemented below. The separate preparation channel is also tested
through the protected API and real TLS with a finite test-only Julia child.
The durable job service and protected submission/result API are connected;
neither transport fixture certifies effective numerical container/native quotas.

### Managed agent service

From `playground/`, validate the operator-owned agent configuration and render its
service definition without installing or starting it:

```sh
./lcm runtime check-agent --config /absolute/path/agent.toml
./lcm runtime agent-unit --config /absolute/path/agent.toml
```

Install the rendered text as the user service `lcm-agent-WORKER_ID.service`, where
`WORKER_ID` is the configured worker identity. Installation, user-manager lifetime
and any required host controller delegation are operator actions. After installing
the file, the operator uses `systemctl --user daemon-reload` and
`systemctl --user start lcm-agent-WORKER_ID.service`. Stop that same service with
`systemctl --user stop lcm-agent-WORKER_ID.service`. No persistent service is
installed or enabled by the LCM CLI itself.

The generated `ExecStart` invokes `lcm runtime start-agent` through the exact Julia
binary/project. Running this action directly from a terminal is rejected: before
opening resource ownership, the agent checks the real user-systemd unit, main PID,
invocation identity, cgroup, fixed commands and shutdown policy. Environment flags
are not accepted as a substitute. The check requires the typed properties used by
systemd 252; unsupported managers fail closed. Regenerate the unit after moving
the checkout/configuration or changing its Julia installation.

`ManagedScientificDriver` inspects the locally installed subset of approved
scientific native and container profiles. Missing environments/engines or mandatory limits are
reported as unavailable; images are never pulled and no executor is prepared on
startup. An empty eligible subset can still run the control scheduler. Native-only
agents do not require or probe a container engine. Both backends share one resource
journal, capacity accounting and lease-owned handles; `ContainerScientificDriver`
is retained as an alias, not a second owner. `ManagedResourceDriver` names the
same physical owner. The CLI composes it through `ManagedAgentResources`, whose
scientific and terminal partitions borrow one journal and one agent lease ledger.
Terminal profiles use the same verified managed-service, host-limit and installed
image admission checks. A blanket `terminal_acceptance_pending` flag no longer
hides provisioned profiles. Eligibility is not readiness: a live owner-fenced
lease and the actual guarded REPL startup marker are still required. Terminals
never fall back to a trusted native process.

The root observes the exact verified systemd incarnation's stop state as well as
SIGINT. This allows cooperative cleanup even if another asynchronous task catches
the signal. It does not treat a replaced or missing service as continued authority.
The finite forced-stop policy still covers the complete service cgroup. Its single
`ExecStopPost` invokes:

```sh
./lcm runtime recover-agent --journal /absolute/original/scratch/resources --worker WORKER_ID
```

These arguments are captured when rendering the unit, not reread from a mutable
configuration during cleanup. Recovery uses no broker credential and restores no
lease or preparation. It does not create an absent journal. A running agent's
kernel journal lock prevents competing recovery. A stopped owner's receipts are
removed only after the exact engine scope/resource identity and absence have been
checked; unresolved resources retain their receipts. Native-unit receipts use the
same dispatcher with the invocation and kernel-group checks described below.
Startup recovery must finish before
any announcement. Automatic restart is disabled; restart is an operator action.

The opt-in real service test creates one unique transient **control-only** unit,
verifies its service identity and exclusive journal, force-kills that exact main
process, checks successful post-stop recovery and removes its failed unit record:

```sh
julia --startup-file=no --compiled-modules=existing --threads=2 \
  --project=runtime runtime/test/managed_agent_systemd.jl
```

Optionally set `LCM_TEST_STOPPED_CONTAINER_IMAGE` to an already cached digest-pinned
image. On Podman the test additionally creates one owned **stopped** container,
leaves its pre-bind receipt for actual agent startup recovery, and verifies the
original container inventory is restored. No container process is started and no
image is pulled. Private temporary diagnostics are retained. This verifies real
control-agent death and stopped-resource recovery, not death with running
scientific/terminal containers or effective executor quotas.

### Lease-owned scientific execution

`ScientificResources(driver)` owns the scientific partition. Its constructor
does not launch anything. `AgentService(config, resources)` binds it to that exact
agent ledger; `recover_owned!` must complete before accepting explicit work.
The driver has six required, operator-owned hooks:

- `installed_profiles(driver)` returns locally verified scientific definitions.
- `recover_owned!(driver)` reconciles only its recorded physical resources.
- `verify_executor!(driver, profile, fence)` rechecks source and resource policy.
- `executor_for!(driver, profile, fence)` supplies the exact owned process supervisor.
- `release_owned!(driver, fence)` returns true only after physical cleanup.
- `close(driver)` joins whole-driver teardown, including partial acquisitions.

There is no default driver which treats a source digest or a process handle as an
isolation-preflight pass. Both production backends use fixed launch policies and
mandatory guarded entry, described below. The
test-only native driver lives under `runtime/test/`, is not selected by
the CLI, and makes no memory/CPU/PID/scratch-isolation claim.

The shared owner exposes the narrow internal orchestration API:

| Operation | Contract |
|---|---|
| `prepare_assigned!(resources, fence, inputs)` | Explicit, single-flight representative preparation; returns a task |
| `execute_assigned!(resources, assigned_job)` | Requires the live fence and exact prepared executor ID, generation and model key; its task returns `ScientificOutput` |
| `prepared_execution(resources, fence)` | Reads the local prepared identity without starting a child; remote callers still require fresh inspection |
| `refresh_preparation!(resources, fence)` | Queries the current child without rebuilding, extending TTL or starting a replacement |
| `cancel_assigned!(resources, fence, request_id)` | Cancels only that pending request, never its successor |
| `scientific_status(resources, fence)` | Returns local activity/preparation dimensions and safe correlation fields |

One request may be pending per assignment. Matching preparations share it; other
work receives a finite busy response rather than an unbounded queue. Separate
assignments own separate processes and model caches. The most recent completed
job task is retained for a matching retry; durable redelivery/result-before-ack
still belongs to the broker owner, not a second in-memory job database.

`AgentJobService`, owned by the root agent, now connects the durable v2 consumer
to that scientific owner through a separate NATS connection. Delivered work
checks for an exact stored result first. A matching result is acknowledged
without touching the child; redelivery without a stored result yields an
explicit `execution_uncertain` failure, not an automatic numerical retry.
The service waits up to five seconds behind a read-only preparation inspection,
then rechecks the target atomically. It does not queue behind another calculation
or preparation. At most the agent's configured capacity is in flight.

Progress acknowledgements keep a long calculation's delivery timer alive. A
terminal acknowledgement follows durable result persistence. If persistence
fails, only the already completed outcome is retried, while the same lease
remains usable and no later than five seconds after the scientific deadline.
Stopping the service cancels its exact pending job requests and joins its tasks;
the existing resource owner remains responsible for physical lease cleanup.

Assigned v2 requests and results now require `PreparedExecution` metadata. The
unreleased earlier v2 shape is rejected; legacy v1 is unchanged. Drain older v2
work before upgrading agents. Never purge v1 streams as part of this change.
The actual child reports its registered result-schema version. Provenance uses
the explicitly labelled `environment-sha256:` digest; it does not invent a
scientific package version. The complete inline broker message is bounded to
256 KiB. With private artifacts configured, result values larger than 64 KiB
use the job-scoped storage path below. Without it, an oversized broker result
reports `result_payload_limit`. The registered scientific consumers use the shared
execution/view component. Service availability is reported by `assigned_execution`;
the complete scientific and deployment acceptance record remains separate.

Lease authority is checked before work, during bootstrap/execution, and before
accepting results. Source/resource preflight is checked before and after work.
The child independently checks preparation/model retention before an operation;
the owner rechecks it after successful execution. Status polling does not extend
model lifetime, and process loss cannot trigger an implicit cold replacement.
Cancellation/lost authority invalidates preparation and retires owned resources.
Controlled operation rejection can retain still-live preparation. Unresolved
cleanup keeps the original occupied handle and requires retry, not early release.

Snapshots distinguish job activity from preparation. Fresh remote reports must
use the child inspection path; a cached local snapshot cannot itself prove an
evictable model still exists. Raw engine log text, private exception contexts and
source paths are excluded from the public status task: only bounded progress,
output-line counts, fixed failure codes and owned identities are exposed here.
The preparation endpoint and shared browser client expose those bounded status
fields through the separate scientific channel. The shared diagnostics also
correlate safe job IDs, executor generations and fixed job-stage codes.

### Private scientific results

Both the coordinator's **control file** and the host's **agent file** accept
the same optional `[artifacts]` table. Use a dedicated private bucket/prefix;
never reuse a public artifact prefix or upload directory. Merely checking a
configuration performs no storage request.

For remote agents, provision an S3-compatible store and use matching locations
but **different** credential files on the coordinator and each worker:

```toml
[artifacts]
backend = "s3"
endpoint = "https://objects.example.org:9000"
bucket = "lcm-private"
prefix = "runtime-v1"
region = "us-east-1"
credentials_file = "secrets/artifacts.toml"
# ca_file = "certificates/ca.pem" # optional private-network CA
```

Each mode-0600 credential file contains only `access_key_id` and
`secret_access_key`. The endpoint accepts a scheme, host and optional port, not
credentials, paths or query strings. Remote transport requires verified HTTPS.
An explicit `allow_loopback_plaintext = true` exception is limited to literal
loopback development/test addresses. No proxy environment, redirects, cookies
or hidden transport retries are used. The existing AWS SDK signs requests.

Provision these distinct object permissions; configuration does not install IAM:

- Coordinator: `GetObject` under `runtime-v1/` and bucket-list permission needed
  to distinguish an absent object from access denied. No write/delete permission.
- Worker `worker-a`: `PutObject` and `DeleteObject` only under
  `runtime-v1/workers/worker-a/`. No result read or other worker's prefix access.
- Anonymous/browser users: no bucket access. Scientific child processes receive
  neither storage credentials nor the host agent's configuration.

Objects use
`PREFIX/workers/WORKER/runs/RUN/BOOT/LEASE/GENERATION/JOB/{sha256,metadata}/DIGEST`
(metadata adds `.json`). This retains the shared content-addressed format but
adds the complete job/assignment scope. Knowing a digest cannot select another
job's objects. The result's legacy-shaped `retrieval_reference` is passive wire
metadata only; the private gateway and browser client never follow that URL.

For trusted single-machine operation, both owners may instead use the **same**
dedicated private filesystem root:

```toml
[artifacts]
backend = "filesystem"
root = "/absolute/private/lcm-runtime-results"
```

Directories are mode 0700 and files mode 0600; symlinks, foreign ownership and
linked/replaced files fail closed. Filesystem mode requires the same trusted OS
identity and filesystem access. It is not a remote-worker substitute for S3.

The host agent serializes successful JSON values and uploads values above
64 KiB; metadata is committed last. Storage accepts at most 4 MiB, while the
current scientific executor's framed IPC retains its tighter **1 MiB complete
response** bound. This change does not promise unlimited scientific output.
An upload failure produces `artifact_unavailable`, not success or an automatic
rerun. A failed remote write attempts finite deletion of its two exact objects.
Filesystem writes synchronize a private temporary file and atomically rename it;
normal failure removes that staging file. The agent holds no disk upload staging
when using S3.

Successful results are intentionally persistent, not deleted on a browser
disconnect or lease release. Provision private-prefix retention for orphaned
objects and a retention period compatible with job history. If remote deletion
is unavailable, cleanup cannot be guaranteed immediately; expiry is the
backstop, not a false successful-cleanup report. The isolated test store uses
two-day expiry; production retention is an operator decision. Filesystem crash
residuals likewise require a stopped-owner maintenance policy, never blanket
deletion of an active artifact root.

Downloads require the authenticated owner, a saved receipt and matching durable
result. The gateway then checks metadata, exact length and SHA-256 before serving
JSON with `no-store`, `nosniff` and attachment disposition. Cross-owner requests
fail before contacting the broker/store. At most four artifact reads run per
coordinator, separately from inventory and status; shutdown joins them. Each S3
request has a five-second deadline and a capped response buffer. The shared
browser `jobArtifact(job_id)` allows a separate 15-second read and never follows
a worker URL. A missing/expired object is explicitly unavailable, not silently
recomputed.

The protected job API retains the original UUID, inputs, deadline and target on
an explicit identical retry. Uncertain jobs cannot be resubmitted by a GET or
reconciliation. Cancellation is persisted before delivery where necessary and
only a durable canceled result establishes its outcome. SQLite **schema 5** is
required; migration remains explicit with the coordinator stopped and backed up.

### Container host prerequisites

The read-only operator check is available before configuring or launching an agent:

```sh
./lcm runtime check-host
./lcm runtime check-host --runtime podman
./lcm runtime check-host --runtime docker
```

Engine discovery is shared with the existing deployment CLI in
`common/container_engine.jl`. Auto prefers actual Docker Engine, recognizes a
Docker command implemented by Podman, and then selects native Podman. An explicit
choice never silently switches engines. The agent path does not require Compose.

The check uses bounded CLI processes: default capacity eight, ten-second command
deadlines and one MiB combined stdout/stderr per command. It never invokes a shell
or inherits stdin. Cancellation, output overflow and timeout retire the original
process and join its readers. Unresolved cleanup remains occupied and can be
retried; private command output, arguments and exception contexts are not status
messages. These are host-CLI bounds, not container resource limits.

Only local Linux engines with the reported cgroup-v2 CPU, memory and PID
prerequisites and seccomp support pass. Docker also requires swap-limit support
and an inspected Unix-socket endpoint; lifecycle commands can use that explicit
socket without a later context change redirecting them. The CLI environment
excludes broker/storage/proxy credentials and remote-engine overrides. Operator
home/runtime-directory references exist only for the host CLI's configuration,
not as mounts or environment for user code.

The current development host returns exit **2** with
`cpu_controller_missing`. This is an expected unavailable result, not a reason
to omit a required limit. Podman documents the dependency of rootless CPU quotas
on host/controller support in its
[run reference](https://docs.podman.io/en/latest/markdown/podman-run.1.html).
Docker's prerequisite fields follow the
[Engine info schema](https://docs.docker.com/reference/api/engine/version/v1.52/).

Passing this check is **not** an executor attestation. Production physical drivers
must still validate the approved image/source, create an exactly owned resource,
verify its effective kernel limits and writable mounts, and recover/retire it
before reporting preparation or release. This command does not pull images,
create containers, certify scratch quotas or confer execution authority.

### Physical acquisition receipts and recovery

`ResourceJournal` stores bounded private JSON records in a dedicated mode-0700
directory. Use persistent local storage for a deployed agent, not an ephemeral
temporary directory or a shared upload location. A kernel-held exclusive lock
prevents two agents from owning the same journal. Files must be regular,
single-linked, mode 0600 and owned by the effective user. Linked paths, replaced
locks/directories, malformed records, FIFOs and unknown entries fail closed.

The physical driver commits `reserve_resource!` **before** its launch command,
then uses `bind_resource!` to record the full container ID or native systemd
invocation ID. The immutable receipt includes the exact assignment fence,
generated UUID name and inspected supervisor/engine scope. Bindings cannot switch
to another physical resource. These records are cleanup identity, not another
lease database and never evidence of retained preparation.

Writes use a private temporary file, file synchronization, same-directory atomic
rename and directory synchronization. These steps follow the
[Linux fsync contract](https://man7.org/linux/man-pages/man2/fsync.2.html).
Startup may remove private incomplete-write fragments only after every committed
record and filename validates. A missing/corrupt ownership marker or foreign
entry requires operator attention; it is not permission to erase the directory.
Closing a journal releases its lock but preserves its records for recovery.

`container_scope` inspects the pinned local engine command and combines its
identity with machine ID and effective UID. Docker uses its daemon ID; Podman
uses its persistent graph-root identity. A changed engine, storage root or context
cannot prove that an old resource disappeared. Restore the original engine
configuration to reconcile that receipt; do not delete it to make startup pass.

`remove_owned_container!` requires that scope, the generated name, every owned
label and any already-bound full ID to match. An interrupted pre-bind acquisition
is resolved by name and labels, then bound to its actual full ID before removal.
Stop/removal targets only that ID. A failed inspect or successful remove command
does not establish absence: a working engine's bounded inventory must confirm
absence in the same scope before the receipt is forgotten. Cleanup remains
available if launch prerequisites such as CPU delegation have since disappeared.

`recover_containers!` attempts all records for its engine even when one fails,
and reports unresolved recovery without touching mismatched resources. Other
backends' receipts remain for their own adapters. The root physical driver must
reconcile **every** backend before announcing availability; this helper's empty
subset is not a whole-agent recovery claim.

Read-only real-engine verification:

```sh
LCM_TEST_CONTAINER_RUNTIME=podman julia --startup-file=no --project=runtime runtime/test/container_scope_local.jl
```

The stopped-acquisition gate requires an explicitly supplied, already cached,
digest-pinned image; no pull is performed:

```sh
LCM_TEST_CONTAINER_RUNTIME=podman \
LCM_TEST_STOPPED_CONTAINER_IMAGE='registry/image@sha256:REPLACE_WITH_CACHED_DIGEST' \
julia --startup-file=no --project=runtime runtime/test/container_stopped_recovery.jl
```

This last test creates a uniquely named **stopped** container with CPU, memory and
PID limits configured, leaves its receipt at the pre-bind crash gap, reopens the
journal and removes only that verified ID. It never starts a container process.
It can therefore test partial-acquisition cleanup on this CPU-constrained host,
but it cannot certify running-worker limits. On unresolved failure it preserves
the private journal/diagnostics instead of discarding the cleanup target.

The journal, container/native recovery, managed agent and shared scientific-driver
integration are implemented. Actual successful quota-enforced native/container
launch verification, remote scientific dispatch and terminal lifecycle integration
remain required. The current host does not delegate CPU to the user manager.

### Native service ownership and recovery

`native_scope(runner)` verifies the fixed local Unix user-bus socket, its owner,
the bus machine identity and the service manager's effective user. It hashes the
machine/user namespace, not a mutable context name or inherited remote bus address.
The durable receipt binds each resource to its generated `lcm-exec-UUID.service`
and exact systemd invocation. `native_resource_description(receipt)` derives its
fixed ownership marker from the shared receipt labels; it grants no lease.

`inspect_native_unit` requires matching transient metadata, invocation and fixed
shutdown policy under the user's `app.slice`. It never loads or starts a unit.
`remove_owned_native!` binds a started pre-bind crash-gap receipt before requesting
stop. A verified inactive, never-started unit has no invocation yet: its exact
owned queued start can be canceled while the receipt remains unbound. A bound
receipt never reverts to this state. Recovery rechecks scope/identity before
mutations and never signals a saved PID. It may
reset only the exact completed owned failed-unit record, not all failed services.

Successful cleanup requires no pending unit job **and** an empty or removed exact
kernel cgroup on an actual cgroup-v2 filesystem. The unit must be absent or have
receipt-matched inactive state and no process. Its inactive definition may remain
cached while another unit references it; recovery never removes that other unit
to force garbage collection. A successful
stop command, failed inspection or missing manager is not enough. A replacement
invocation, queued start, surviving descendant, changed scope or incompatible
cleanup policy retains the receipt. Polling and each manager command are bounded;
an unresolved teardown may be retried without replacing its ownership identity.
`recover_native!` attempts all native receipts and leaves container receipts to
their own adapters. The root `recover-agent` dispatcher reconciles both kinds.

These cleanup functions do not grant scientific admission or attest limits.
The separate native launch policy below uses the same receipts and enforces its
own preflight and child-entry checks. Arbitrary Julia terminal input remains
container-only.

From `playground/`:

```sh
julia --startup-file=no --project=runtime runtime/test/native_recovery.jl
julia --startup-file=no --compiled-modules=existing --threads=2 \
  --project=runtime runtime/test/native_recovery_systemd.jl
```

The opt-in second gate creates only receipt-owned `sleep` services, with an
automatic 120-second maximum runtime. It removes one while verifying another is
unchanged, exercises the actual agent recovery dispatcher after journal closure,
and cancels a queued dependency before its process starts. The owned start barrier
is finite and also cleaned up. No numerical process, terminal or image is started. The test
preserves its private journal/diagnostics if cleanup cannot finish. It does not
relax any executor admission requirement to test the cleanup mechanism.

Verification (from `playground/`):

```sh
julia --startup-file=no --project=runtime runtime/test/runtests.jl
julia --startup-file=no --project=runtime runtime/test/scientific_profile_processes.jl
```

The second test uses the actual line and power-flow environments with real elapsed
time and independently maintained lease renewals. It verifies concurrent numerical
work, preparation reuse and independent release through the shared owner. Its
explicit native test driver is not evidence of a production physical-resource
preflight. The service-level `assigned_execution` flag is not an attestation of
this test driver's resource isolation or a substitute for production preflight.

### Trusted native scientific launch

From `playground/`, inspect native prerequisites without creating a resource:

```sh
./lcm runtime check-host --runtime native
```

The local user-systemd backend requires non-root identity, cgroup v2 with delegated
CPU, memory and task controllers, and user namespaces. The CLI reports fixed
failure reasons and exits nonzero when prerequisites are missing. Passing this
read-only check does not establish effective limits or readiness. On this host it
exits 2 with `cpu_controller_missing`; no alternative quota or unconfined fallback
is substituted.

`NativePolicy` accepts only approved `scientific`/`trusted_process` profiles.
The shared driver rechecks its managed-agent incarnation, local scope and source
fingerprint, records an acquisition intent, then constructs the fixed
`systemd-run --user --pipe --wait` command. Its generated service is bound to the
agent service's lifetime and uses the recovery policy above. Standard IO carries
the existing bounded scientific frames, not a terminal PTY.

The service requests CPU quota, memory with no additional swap, task count,
read-only ordinary filesystems, private users/devices/network/IPC and finite
descriptor/core/message-queue limits. Writable `/tmp` and `/dev/shm` are private
tmpfs mounts whose combined capacity cannot exceed the approved scratch budget.
CPU and scratch rounding is downward. No persistent scratch directory is created
for an executor; its mounts disappear when the owned service group is retired.

The fixed command uses `env -i`: broker credentials, proxy settings, inherited
Julia load paths and user-bus addresses are not forwarded to the scientific
process. Only explicit offline Julia settings, installed depot paths, guard
expectations and approved package UUID/name are supplied. Source/binary/depot
paths must be absolute and visible outside `/tmp` and `/dev/shm`; paths containing
control characters, `:`, `$` or `%` are rejected rather than expanded by the
service manager. Install dependencies before admission; launch never instantiates
or downloads packages.

`worker/core/src/native-guard.jl` runs before the scientific package import. It
checks actual kernel identity, the exact generated cgroup, finite quotas, zero
swap, mount policy/scratch capacity, no-new-privileges and dropped capabilities.
The kernel checks share `ExecutorLimits` and quota validation with the container
guard (`ContainerLimits` remains a compatibility alias). Failed guard entry exits
78 with a fixed short diagnostic and does not enter the profile. These checks are
required because [systemd documents that some sandbox settings can be unavailable
or silently ineffective](https://raw.githubusercontent.com/systemd/systemd/v252/man/systemd.exec.xml).

After successful preparation/execution, the driver binds the exact service
invocation and rechecks its command, agent dependency and configured policy before
accepting the result. A replacement process cannot inherit readiness. Cancellation
or failed launch retains the receipt until the original attached process and the
exact service/group are both retired. Native profiles remain **trusted code**:
read-only host visibility is not an arbitrary-code security sandbox.

Current evidence covers policy/drift tests, real denial before evaluation,
native acquisition without start, and real control-agent crash/recovery while
native prerequisites are missing. A successful numerical launch with effective
limits still needs a correctly delegated host. The opt-in negative admission gate
on a host with unavailable native prerequisites is:

```sh
LCM_TEST_NATIVE_UNAVAILABLE=1 julia --startup-file=no --compiled-modules=existing \
  --project=runtime runtime/test/managed_agent_systemd.jl
```

It launches only uniquely owned control-agent test services, then forces their
exit and verifies post-stop recovery. It does not start a scientific or terminal
executor, install a persistent service, or relax host admission requirements.
For inert validation of the generated native service directives, run
`julia --startup-file=no --project=runtime runtime/test/native_policy_unit.jl`.
This uses `systemd-analyze --user verify` on private generated files, not a
successful service launch or effective-quota certification.

## Shared container policy and image entry

`ContainerPolicy(profile, receipt)` is the single source for scientific and
terminal container configuration. `create_owned_container!` checks the selected
host, pinned local image and exact receipt scope, creates only a stopped resource,
then checks its entry command, environment, ownership, namespaces, privileges,
mounts and resource limits. It never pulls, starts or automatically retries an
image. Failed acquisition retains a recoverable intent; its caller must complete
physical cleanup before releasing the assignment.

The immutable image entry loads the lightweight `verify_container_isolation`
check before numerical imports or the Julia REPL. Kernel evidence and configured
limits are distinct from preparation and lease authority. Missing CPU delegation
remains an error, not an affinity-only fallback. The shared policy explicitly
accounts for Docker/Podman schema differences; individual profiles do not provide
their own flags or security rules.

See [approved executor images](../worker/containers/README.md) for the three
separate build targets, source-layout contract, operator provisioning boundary,
mandatory limits and current verification status. No broker or host credentials
are copied into those images. The managed scientific driver rechecks the actual
service, engine, image and container policy for its exact assignment, retains
partial acquisition, and owns the attached CLI plus physical removal. An attached
CLI cleanup failure does not suppress the physical-removal attempt; unresolved
process/receipt ownership keeps capacity occupied. Remote scientific dispatch is
connected through the protected job API. The shared scientific controls are
available, but their real deck/workbench consumers and the terminal relay remain
incomplete, so their public capabilities stay disabled.

From `playground/`:

```sh
julia --startup-file=no --project=runtime runtime/test/container_policy.jl
julia --startup-file=no --project=runtime runtime/test/container_image_layout.jl
julia --startup-file=no --project=worker/core worker/core/test/container_isolation.jl
julia --startup-file=no --project=worker/core worker/core/test/container_entry.jl
```

The optional `runtime/test/container_stopped_policy.jl` requires
`LCM_TEST_STOPPED_CONTAINER_IMAGE` to name an already cached digest-pinned image.
It creates only stopped acquisition metadata, checks both scientific/TTY command
policies and removes only its new receipt-owned resources. It deliberately does
not call production image admission or start the base test image. This is a real
engine configuration/recovery gate, not a running-container isolation result.

## Terminal PTY ownership — internal transport

`src/TerminalProcess.jl` provides the bounded local transport used by the managed
private-terminal driver. `TerminalProcess(; limits=TerminalIOLimits())` is passive;
the physical supervisor must retain it before `start_terminal!` receives its fixed,
operator-built container-attachment command. The command needs an explicit
environment. It is not accepted from a browser, and this internal transport is
neither a new native-terminal CLI mode nor an isolation substitute.

The implementation uses a nonblocking Linux PTY master and the FileWatching
standard library. A TTY-enabled Docker attachment requires an actual terminal,
not redirected stdin. The pump avoids `uv_tty_t` master writes because libuv may
fall back to blocking writes for a PTY master. See the
[Docker attachment contract](https://docs.docker.com/reference/cli/docker/container/attach/)
and [libuv TTY implementation](https://github.com/libuv/libuv/blob/v1.x/src/unix/tty.c).

Current transport guarantees:

- One bounded I/O pump per handle, a fixed output ring, a maximum input queue and
  chunk size, and output-rate, stalled-write and command-lifetime deadlines.
  Continuous output yields to other scheduler tasks.
- `read_terminal` returns raw bytes, a byte cursor, the latest sequence and an
  explicit gap when older bytes were evicted. UTF-8 and escape sequences may span
  chunks; the transport neither interprets nor executes them.
- `write_terminal!` acknowledges local queue admission only. It does not claim
  Julia evaluated input. Sole-writer authority, stream sequence checks and
  reconnect policy belong to `TerminalResources`, described below.
- Resize uses the PTY dimensions and signals only the original attached command.
  The container CLI must forward terminal behavior to its own PTY; this
  local test alone does not certify Docker/Podman resize propagation.
- Exit, overflow, stalled writes and explicit close retire the original command
  handle, join the I/O watcher, then close the descriptor. A stuck task retains
  ownership for a bounded cleanup retry. Descriptors cannot be recycled while a
  live watcher still references them. Attached CLI exit is not container removal.
- Input chunks are not durable or logged. Ordinary displays of process/output
  handles redact their private bytes and command details.

Repeat the focused test from `playground/`:

```sh
julia --startup-file=no --threads=2 --project=runtime runtime/test/terminal_process.jl
```

This launches **finite trusted local fixtures**, including an actual Julia REPL
with a separate kernel deadline. It verifies multiline input, history, completion,
Unicode, resize, interrupt, independent REPL state, floods, stalled input,
cleanup retry and repeated descriptor retirement. It starts no container and
does not test effective CPU/memory/PID/scratch isolation. Public `private_terminal`
reports whether the configured relay is available, not whether an eligible worker
or ready REPL exists. Missing delegated CPU still fails closed for real terminal
containers. The opt-in `runtime/test/physical_terminal.jl` and
`runtime/test/run-broker.sh physical` exercise actual resources and the managed
path separately; see the progress ledger for their results and outstanding gates.

### Shared physical owner and private terminal sessions

`ManagedAgentResources` owns one physical `ManagedResourceDriver`, a
`ScientificResources(ManagedScientificView(parent))` partition and a
`TerminalResources(ManagedTerminalView(parent))` partition. The existing
`AgentService` remains the sole control scheduler and lease ledger owner. Its
scientific and job services receive only the scientific partition through
dispatch. There is no competing journal or separate terminal capacity pool.

Borrowed partitions reject cross-kind acquisition/release. Closure first stops
admission and joins in-flight acquisitions, then retires that partition's exact
handles. It does not close its sibling or the shared runner/journal. Root closure
revokes both partitions and attempts both cleanups and physical-parent closure;
unresolved ownership remains retryable. Failed acquisition occupies the same
total capacity until its cleanup succeeds.

`AbstractTerminalDriver` requires verified profiles, recovery, acquisition,
fixed startup, acquisition-bound marker retrieval, physical verification,
exact release and close. The production adapter creates only a stopped,
policy-verified container, retains its receipt and PTY, and attaches using its
full immutable container ID. The successful entry guard installs a fixed Julia
REPL initialization hook with a receipt-bound marker. `TerminalResources` waits
for that marker, rechecks physical policy and process liveness, and removes the
startup prefix from exposed output before accepting writer input. This is startup
evidence, not a claim that later Julia commands will succeed or finish.

The session API is internal to the authorized runtime path, not a public endpoint:

- `open_terminal!` requires a usable assignment and a writer UUID. Matching
  retries return the same stream; another writer is rejected. Ended streams do
  not implicitly restart. Replacement must be an explicit higher-level action.
- `terminal_status` exposes bounded state and cursor counters, never the writer
  token, input, output or command. `read_terminal` checks exact run/lease/stream
  identity and returns bounded raw bytes with relative cursors and gap flags.
- `write_terminal!` checks connected-writer authority and consecutive sequence
  numbers. Only an exact retry of the latest admitted chunk is acknowledged
  without replay. An altered, older or skipped sequence fails. A queue ACK is
  not an evaluation ACK. No raw terminal input enters durable scientific jobs.
- `resize_terminal!` changes cell dimensions. Ctrl-C uses sequenced input byte
  `0x03`, preserving ordering rather than creating another input path.
- `keepalive_terminal!` refreshes connected-writer presence, not idle time or
  lease authority. The default presence interval is 15 seconds. Its expiry starts
  disconnect grace at the original deadline even when worker control remains
  healthy. Cached heartbeat replies do not refresh presence; only explicit open
  by the same writer may reconnect within grace.
- `disconnect_terminal!` begins grace once. Only the same writer can reconnect;
  repeated disconnects and passive output/status reads do not renew it or idle
  time. Defaults are 120 seconds startup, 30 seconds disconnect grace and 1,800
  seconds idle, with a 15-second task-join cleanup attempt. PTY lifetime, output
  rate and byte bounds remain independently enforced by `TerminalIOLimits`.
- Lease loss, session deadlines, failed startup and natural exit stop transport
  and attempt exact physical retirement. A temporarily refused container removal
  does not leave the byte stream running. Failed cleanup retains ownership; ended
  session records remain until lease release, not as remembered live REPL state.
- Explicit stop retains that ended identity. Restart joins exact cleanup before
  opening a fresh stream/Julia namespace. A failed cleanup cannot allocate a
  replacement. The relay preserves the original request revision on an exact
  restart retry, so an uncertain reply cannot cause two replacements.

The regression tests use the real Julia `--interactive` startup and REPL under a
finite, explicitly test-only native driver. They verify startup-marker gating,
separate Julia state, writer exclusion, non-replayed input, interruption,
reconnect, deadlines, lease loss and retryable cleanup. Production advertisement
includes only installed terminal profiles admitted by the managed host's actual
image/isolation checks. The separate physical engine gates verify those checks
and fault containment; no browser/CLI native-terminal fallback exists.

### Private terminal relay and browser socket

The protocol's `TerminalCommand` and `TerminalReport` are transient version-2
records, separate from scientific jobs. They bound raw byte chunks to 8 KiB and
encoded frames to 64 KiB. Each command carries the full assignment fence,
monotonic revision, request UUID, exact stream UUID and, for writer actions,
the sole writer UUID. Input additionally has its own consecutive byte-chunk
sequence. The agent keeps only its latest request digest/report and one finite
action task per occupied lease; completed tasks no longer retain input bytes.

`BrokerTerminal`, `AgentTerminalService` and `TerminalCoordinator` use independent
connections/tasks. Exact subjects include worker identity, worker boot, lease UUID
and generation. Provisioned broker credentials enforce worker and direction;
full run/owner/profile/stream authority is checked again from the payload and
live lease. These are not per-user broker credentials: users never receive NATS
credentials, and the trusted gateway enforces private owner access. No terminal
subject uses JetStream. Lost connections retire their old queues.

The private gateway route is:

```text
GET /runtime/api/assignments/<lease UUID>/terminal  (WebSocket upgrade)
```

Existing origin/proxy-principal checks run before upgrade. Only the actual owner
may attach, including when another principal is an inventory administrator.
There is one live browser attachment per lease. Every received action and every
outgoing report rechecks the full live assignment; stopping a run closes its
socket. A browser cannot choose a fence, revision, environment, container command
or worker by supplying JSON fields.

The socket sends `hello`; its first private input must be
`{"action":"attach","writer_id":"<UUID>"}`. Subsequent action records have
exactly `action`, `request_id`, `session_id`, `input_sequence`, `after`, `columns`,
`rows`, `bytes` and `retry`. Unused values are zero/null/empty under the shared
protocol grammar. Open, status/read, input, resize, keepalive, disconnect, stop
and restart use the same agent resource owner. A returned report is an action
acknowledgement, not an indefinitely valid readiness grant. Read replies carry
bounded raw output and explicit cursor/gap information; consumers request the
next chunk only after consuming the previous one.

Only one request is in flight and one may wait in the browser relay. The receive
codec queue stores four frames, plus the reader's one blocked frame; fragmentation
and compression are disabled. Input admission is rate-bounded. Slow writes expire
after three seconds, initial attachment after five seconds and silent sockets
after 45 seconds. Consumers send fresh keepalive actions every five seconds,
including during REPL startup. Browser close sends a bounded best-effort
disconnect; worker-side presence/lease deadlines cover lost cleanup messages.

Intentional **Disconnect** is serialized after any current request and waits for
the remote writer-release acknowledgement. Input and reconnect remain disabled
while it is pending. A confirmed disconnect is not sent a second time during
socket cleanup. Reconnect may wait at most five seconds for an already-closing
predecessor's exact attachment slot; active writers are still rejected
immediately, and ownership/lease checks continue during the wait. Expired cleanup
does not free or steal the old slot and does not trigger an automatic retry.

`TerminalSockets.jl` contains a narrow HTTP **2.6.6** compatibility adapter because
that release otherwise creates an unbounded receive channel. It installs a bounded
channel before the reader starts, without changing any global HTTP method or
package. A different HTTP version or already-started reader fails closed. The
adapter's flood test verifies queue bounds and joins a blocked reader on teardown;
HTTP upgrades require reviewing this adapter and rerunning that test.

Input, restart, resize, stop and disconnect are never automatically replayed.
An `uncertain` response retains the latest request identity and digest; an explicit
exact retry resends its original command revision. A later query may supersede
an expired request, after which the old retry is rejected. The browser must not
automatically retype uncertain input. Agent-side input sequencing still guards
against repeated queue admission. Terminal bytes remain outside SQLite, ordinary
diagnostics, scientific results and X-ray.

From `playground/`, the isolated TLS integration harness is:

```sh
CONTAINER_RUNTIME=podman bash runtime/test/run-broker.sh terminal
```

It uses its own temporary NATS container and finite, explicitly trusted native
REPL fixtures; it removes only those owned resources. It is not a verification
of numerical/terminal container isolation. Actual Docker and quota-enforced
Podman launch are exercised by the separate [physical gate](PHYSICAL_ACCEPTANCE.md)
on a capable host; the ledger records the completed owned-Kubuntu matrices.

### Reusable Julia terminal view

The same component is rendered by Bonito in gallery, workbench and live-frame
contexts; a plain owned browser page can use its shared assets directly. There
is no page-specific terminal implementation or separate theme palette:

```julia
using Bonito, UUIDs, LineCableModelsPlayground

client = RuntimeClient(UUID(ENV["LCM_RUN_ID"]))
terminal_panel = DOM.div(
    WorkerSelector(client, :terminal; profiles=("julia-terminal",)),
    JuliaTerminal(client, :terminal; title="Julia REPL", rows=18),
)
```

Return `terminal_panel` from an owned branded live view or compose it into the
workbench workspace. The outer shell supplies the shared brand/control styles;
the component does not install a second global palette.

The application definition must declare the `terminal` role and the operator
must approve its container-backed profile. Passing a role or run UUID grants no
access. The shared client does not poll terminal roles for scientific preparation,
and the general control panel suppresses scientific preparation actions for those
roles; terminal startup belongs to Connect. `/widgets/julia-terminal` is an
inventory-only public specimen, so its Connect action is disabled. A configured
private relay is a capability, not readiness or profile eligibility: Connect
still requires the current owner, a usable assignment, an approved immutable
image and successful host/container admission. The component cannot bypass
those checks. Physical acceptance evidence is recorded separately in the ledger.

`src/widgets/JuliaTerminal.jl` contains the passive Julia value and safe X-ray
descriptor. `assets/runtime-terminal-client.js` owns the serialized private
socket and `assets/runtime-terminal.js` owns xterm rendering and teardown.
`runtime-terminal.css` owns only layout; shared brand/control tokens govern both
themes. Pinned vendor assets are local, with their build and licenses recorded
in `assets/vendor/README.md`. Raw input/output never becomes a Bonito Observable.

Rendering and slide entry connect to no REPL. Connect is explicit and reports
starting before input becomes available. The client sends fresh keepalives,
coalesces resized cell dimensions, and permits at most one socket request at a
time. Its input queue is at most 32 KiB, each private chunk is at most 8 KiB,
and an oversized paste is rejected whole. No output poll occurs while xterm is
still consuming the previous chunk; a stalled render has a three-second bound.
Retained output gaps reset the display parser before showing the available tail.
The visible scrollback is limited to 1,000 lines.

Interrupt sends Ctrl-C and discards unsent input. Stop and restart ask for
confirmation; restart waits for worker cleanup before obtaining a new Julia
namespace. Clear view affects local scrollback, not Julia variables. Intentional
disconnect discards unsent input and waits for writer release before closing.
Page departure, unmount, changed assignment or stale inventory close the socket
immediately and discard unsent input. The worker owns bounded grace/idle/lease
cleanup even when no disconnect acknowledgement can reach the browser.

A lost reply never causes input or mutations to be replayed. Reconnect inspects
the surviving stream and its admitted input sequence. After an uncertain action,
output remains readable but **Resume input** is explicit; inspect it or choose
a clean restart before re-entering code. An input acknowledgement means queued
bytes, not completed evaluation. The writer identity exists only in the mounted
component's memory: page reload cannot claim restoration of a lost session.

Terminal keyboard events do not bubble into presentation navigation. User text
selection/copy remains available. Output-controlled OSC title, hyperlink,
clipboard and palette changes are consumed without applying them, window
operations are disabled, and output is never inserted as application HTML.
Ordinary Julia ANSI/VT output remains terminal data. X-ray exposes only role,
title, row configuration and action names—not bytes, writer identity or history.

Focused browser/component checks from `playground/`:

```sh
node test/integration/runtime_terminal.mjs
node test/integration/runtime_terminal_browser.mjs
bash runtime/test/run-bonito.sh
```

The first uses explicit socket fixtures; the second uses real Chrome/xterm with
explicit HTTP/socket fixtures. The Bonito gate verifies actual proxied asset
mounting, theme parity and X-ray in separate owned UI hosts. The terminal TLS
suite additionally exercises a real leased native REPL fixture; none of these
native fixtures asserts production container isolation.

## Protected scientific job API

`ControlService` owns one `JobCoordinator`, sharing the existing lease,
preparation and SQLite authorities. Its job/result connection and bounded tasks
are separate from heartbeat, preparation and HTTP handling. Construction is
inert; starting control starts reconciliation, not a scientific computation.

All endpoints below use the existing authenticated same-origin gateway. POST
requests additionally require its Origin and `X-LCM-Request: 1` checks. Missing
and foreign resource IDs produce the same 404; ownership precedes parsing input
or accessing the broker. The browser cannot supply a worker, executor, image,
command, lease generation, job ID or deadline for execution.

| Method and path under `/runtime/api` | Meaning |
|---|---|
| `POST /assignments/<lease>/jobs` | Persist a registered operation, passive input object and caller `request_id`; return 202 with the server-authored receipt |
| `GET /runs/<run>/jobs` | At most 256 recent owned receipts, without result payloads |
| `GET /jobs/<job>` | Owned submission, provenance, state and cancellation flags |
| `GET /jobs/<job>/result` | Exact durable result, or `result: null` only when the broker confirms absence; transport failure is 503 |
| `GET` or `HEAD /jobs/<job>/artifact` | Owned, hash-verified private JSON bytes; missing inline-only result is 404, unavailable storage is 503 |
| `POST /jobs/<job>/cancel` | Persist cancellation intent with a caller `request_id`; return 202, not a claim that execution has stopped |

Submission requires fresh explicit preparation and an allowlisted operation on
the acknowledged assignment. An identical retry preserves the original job,
absolute deadline, input hash and prepared executor target, even when that
assignment later expires. Changed inputs with the same request identity fail.
Automatic publication retries retain that immutable request and remain inside a
two-minute window, strictly shorter than the stream's ten-minute deduplication
window. An uncertain computation is never rerun merely to obtain a result.

`queued` and `submitted` describe delivery bookkeeping; `succeeded`, `failed` and
`canceled` require a matching durable result. `revoked` means assignment authority
was lost, not that an interrupt was confirmed. After the job deadline plus a
five-second result grace, an unresolved receipt becomes `uncertain`. A later
owned result read can resolve uncertainty without publishing another job.
Historical results do not establish current readiness or permission to execute.

Cancellation intent and worker acknowledgement are separate from job outcome.
The worker records a bounded tombstone for the exact job and assignment before
interrupting an active operation. A queued canceled job cannot start after
delivery or after the child is re-prepared under the same lease. Its original
delivery still produces an actual canceled result. A completed durable result
wins over a late cancellation. Lost cancellation replies may safely resend that
same tombstone; uncertain preparation requests are still never replayed.

Reconciliation has at most four concurrent finite attempts and rotates pending
receipts fairly. HTTP result reads and artifact reads each have a separate bound
of four. Shutdown joins all three before closing the job connection. Diagnostics use the existing bounded,
owner-filtered event buffer, adding job/executor IDs and fixed stage codes, never
scientific inputs, result values or arbitrary exception text.

The shared browser `RuntimeClient` provides `submitJob`, `listJobs`, `job`,
`jobResult`, `jobArtifact` and `cancelJob`. It validates passive inputs and returned provenance,
preserves explicit retry identity after uncertain/invalid mutation replies, and
never follows a worker-supplied artifact URL. The shared `RuntimeJob` view tracker
and `ScientificJob` Bonito component consume those methods. The registered ICHQP
showcase deck and CableStudy workbench use the same scientific adapters, verified
against separate direct-engine calculations. See [SCIENTIFIC_CONSUMERS.md](../SCIENTIFIC_CONSUMERS.md)
and the [acceptance ledger](../RUNTIME_PLATFORM_PROGRESS.md). Artifact delivery uses the
[private scientific result contract](#private-scientific-results) above.
`assigned_execution` reports the open scientific service; `private_terminal`
retains its independent acceptance gate. Results beyond the configured storage/IPC bounds
still fail explicitly; the legacy public artifact hash route does not serve
private runtime results.

The browser polls the durable job receipt independently of result delivery.
Queued/submitted jobs and revocation remain visible even when no result stream
exists or artifact delivery is unavailable. Successful data is fetched separately,
with bounded requests and retained provenance; a delayed successful result can
recover without replaying its calculation.

On job-connection startup, the coordinator provisions/verifies the bounded
streams for configured scientific workers before reporting that transport online.
It does not allocate a lease, import an engine or prepare a model. First-use
transport compilation therefore precedes live job publication; the publication
path still rechecks stream policy. A conflicting broker policy requires operator
reconciliation, never destructive replacement or relaxed lease deadlines.

## Database migrations

The TOML configuration remains version 1. The SQLite schema is now version 5:
it retains runs, worker registrations/operator actions and assignment fences,
and durable job-submission receipts, and adds owned cancellation intent.
New databases use the current schema; an existing schema-1, schema-2, schema-3
or schema-4 database is never upgraded implicitly
during startup.

A receipt retains the original server-authored job ID, inputs, deadline and
prepared target across retries/restart. It is not proof of successful delivery
or execution; JetStream remains the result authority. Admission permits one
pending receipt and at most 256 receipts over an assignment's lifetime. Reusing
an idempotency key with changed inputs is rejected. History remains owner-checked,
including after a lease is released, and cannot revive that lease. No existing
operator database is migrated automatically.

Stop its runtime, then run:

```sh
./lcm runtime migrate --config runtime/local.example.toml
```

Migration takes the same exclusive lock as the UI supervisor and refuses to run
while that coordinator owns the database. It creates a private
`<database>.schema-<old-version>-backup-*/runtime.sqlite` snapshot beside the original database
before applying the transactional schema change. Repeating the command on a
current database performs no migration and creates no additional backup.

The snapshot uses [SQLite VACUUM INTO](https://www.sqlite.org/lang_vacuum.html),
which copies a consistent database without rewriting the source. Keep the backup
if a migration fails. A valid backup deliberately remains operator-owned recovery
data, not disposable run scratch. Registration history is not evidence of current
worker liveness or a prepared process.

### Backup and restore procedure

Stop admission and the coordinator before a recovery operation. Back up the
SQLite database using SQLite's [`.backup` command](https://www.sqlite.org/cli.html),
not a lone copy of a live database file. Use an existing private backup directory
and a new destination filename. Example paths below are operator-selected, not
default LCM locations:

```sh
umask 077
sqlite3 -readonly /absolute/state/runtime.sqlite ".backup '/absolute/private-backups/runtime-001.sqlite'"
sqlite3 -readonly /absolute/private-backups/runtime-001.sqlite 'PRAGMA quick_check; PRAGMA user_version;'
```

Require `ok` and the expected schema version. Keep the configuration, application
checkout/version, approved profile fingerprints, NATS/JetStream state and private
artifact store backups together in the recovery record. Protect credentials and
encryption keys separately with restricted access. SQLite alone does not contain
artifact bodies, Julia variables or prepared models.

For restore, retain the original database and backup unchanged. In a **new private
directory**, restore into a **new database filename** using SQLite's `.restore`
command, then repeat the checks above. Do not overwrite an active database or
mix another database's WAL/SHM files into the restored copy. Point a copied runtime
configuration at the restored database and run `lcm runtime check`; if its schema
is older, run the explicit migration command before startup. Never run two
coordinators against the same state.

Restart the gateway/coordinator and agents through their normal ownership checks.
Keep each agent's original resource journal: recovery must reconcile its exact
recorded resources before announcing availability. Do not erase receipts to make
startup succeed. Old leases, old UI sessions and a saved `ready` label cannot
restore execution authority; launch a new run, assign and explicitly prepare.

### Operational acceptance and measured preparation

| Process | Julia project / approved environment | Loads scientific engines? |
| --- | --- | --- |
| Gateway, coordinator, host agent | `runtime/Project.toml` | No |
| Bonito UI host | `Project.toml` | No |
| Shared executor protocol/supervision | `worker/core/Project.toml` | No |
| Line-parameter executor | `worker/profiles/line-parameters/Project.toml` | LineCableModels only |
| OHL/UGC executor | `worker/profiles/power-flow/Project.toml` | PowerImpedance and its required numerical dependencies |
| Private terminal | Approved image with `worker/profiles/julia-terminal/Project.toml` | Only what that approved image and its user session load |

Profile approval also requires the declared version, environment fingerprint or
image digest and effective limits. A project directory alone grants no authority.

Use [the aggregate acceptance runner](../VERIFICATION.md#repeatable-commands) on
a dedicated test checkout. Do not edit its scripts during execution or overlap
groups; it serializes aggregate runs with a kernel-held lock. This avoids mixing
publication builds and unrelated compiler-heavy regression workloads into a
timing-sensitive transport measurement. The scenario itself still exercises
separate owners, independent roles, cancellation and worker loss.

On this development machine, the fully passing `Ze6jr9gp` scientific scenario
measured line cold preparation at 24–25 s, corridor cold preparation at 168–179 s,
and repeated explicit preparation at approximately 4.1–4.4 s including status
polling. Explicit cold corridor cancellation/replacement recovery took
approximately 371 s. The complete gate passed 43 browser checks, 32 protected
runtime assertions and 39 separate direct-engine comparisons. These are observed
fixture timings, not startup guarantees or resource-isolation proof.
The test profiles declare 4 GiB budgets per scientific executor, 600 s preparation
and 120 s job deadlines; the finite native fixture does not enforce OS quotas.
Capacity planning for a deployed host still requires measured peak memory and
effective CPU/memory/PID/scratch-limit verification there.

The [implementation ledger](../RUNTIME_PLATFORM_PROGRESS.md) distinguishes each
passed subgate from whole-suite failures and completed host acceptance. The original
IT-managed development machine lacks delegated CPU control and has a Podman
`docker` shim; the separately provisioned owned Kubuntu host passes both physical
engine matrices and separate-computer browser/science/terminal/S3 gates, with
exact teardown verified. These are private rehearsals, not public deployment.
Always run `check-host` on the actual target. Do not remove
mandatory limits or enable an unconfined terminal to make an unsupported host
appear eligible.
