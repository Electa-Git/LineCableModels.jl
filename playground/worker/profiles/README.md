# Scientific executor profiles

These environments separate the v2 numerical dependencies while reusing the
existing operation, cache and executor sources. They are internal child-process
entry points, not public servers or a replacement broker protocol.

| Environment | Numerical dependency | Operations |
|---|---|---|
| `line-parameters/` | Local `LineCableModels` | `cable.geometry_summary`, `cable.constants`, `cable.frequency_sweep`, `line.frequency_scan` |
| `power-flow/` | Pinned `PowerImpedance` revision | `powerflow.prepare`, `impedance.evaluate` |
| `../core/` | None | Shared validation, framing, cancellation, preparation and cache machinery |
| `julia-terminal/` | None | Pinned lightweight environment for the separate container-only Julia REPL |

Neither profile depends on Bonito, NATS, AWS or the other numerical engine. The
legacy `LineCableModelsWorker` remains available for v1 consumers and still has
its original combined registry. Both paths use the same `worker/src/Executor.jl`,
`OperationRegistry.jl`, `Cache.jl`, and scientific adapter files. Keep the complete
worker source tree when packaging a profile; do not distribute only one subfolder.

From `playground/`, install the pinned environments explicitly:

```sh
julia --project=worker/core -e 'using Pkg; Pkg.instantiate()'
julia --project=worker/profiles/line-parameters -e 'using Pkg; Pkg.instantiate()'
julia --project=worker/profiles/power-flow -e 'using Pkg; Pkg.instantiate()'
```

Installation is an operator action, not a page-construction or worker-discovery
side effect. The manifests fix resolved dependencies; local path packages remain
mutable source trees. Agent preflight must verify their complete approved source
fingerprint before advertising them. Merely finding a manifest is insufficient.

Compute the read-only native source digest with the existing CLI:

```sh
./lcm runtime fingerprint --project worker/profiles/line-parameters
./lcm runtime fingerprint --project worker/profiles/power-flow
```

The command prints one SHA-256 value. Copy it into the operator-owned profile's
`fingerprint`; run it on the deployment checkout with the intended Julia version.
It never loads an engine, downloads dependencies, starts a child or contacts NATS.
`verify_native_environment(profile)` repeats the inspection and rejects drift.

This covers the root manifest, all local path dependencies' `Project.toml`, `src/`,
`ext/`, `deps/`, artifact declarations and local preferences. `RuntimeSources.toml`
explicitly covers the shared executor/adapter files outside package directories;
literal cross-package includes are checked against those declarations in tests.
Additional dynamic includes or data outside those trees must be declared there.
Only explicit relative files are accepted, not directory/glob expansion.

The selected `Manifest-vMAJOR.MINOR.toml`, if present, supersedes `Manifest.toml`,
following [Julia's manifest selection](https://pkgdocs.julialang.org/v1/toml-files/).
Alternate `JuliaProject.toml`/`JuliaManifest.toml` naming is currently rejected,
not guessed. The manifest must match the exact running Julia version. Relative
checkouts with identical contents produce the same digest on the same Julia
version/platform; literal absolute paths inside TOML remain part of that digest.

The reader rejects symlinks, missing/unpinned dependencies, source identity/version
mismatches, observed concurrent edits and excessive entries/bytes/depth. It checks
local source contents but trusts the operator's Julia distribution and installed
depot: pinned third-party Git trees identify dependencies, not an integrity audit
of installed artifacts. This is neither a sandbox preflight nor readiness. Never
edit trusted sources under a running prepared executor; retire it and reapprove
the new fingerprint first. Agent launch/recheck integration remains pending.

## Narrow profile contract

`AbstractScientificProfile` requires four dispatch hooks:

- `operation_registry(profile)` returns a passive, closed `OperationRegistry`.
- `validate_preparation(profile, parameters)` validates and normalizes inputs.
- `prepare!(profile, context, parameters)` runs a representative workload and
  returns `PreparedWorkload(evidence; resources=...)`.
- `cleanup!(profile, cache)` releases its retained preparation state.

`PreparedWorkload.resources` identifies model keys required by its evidence. The
common executor rejects reuse if a model expired or was evicted. Line preparation
actually computes the requested geometry at the first and last frequency; it
does not claim to precompile every future Julia specialization. Power-flow
preparation actually constructs, solves and linearizes the existing OHL/UGC case.

The existing `ExecutorSupervisor` accepts an operator-owned profile command.
`prepare_supervised!` uses the same framed channel, cancellation and timeout path
as scientific execution. Bootstrap only means that the command reader is present.
Preparation replies contain normalized input identity, representative-workload
evidence, cache hit/miss/shared status and idle TTL. The owning agent must still
bind that evidence to the live executor, environment fingerprint and exact lease
generation. A replacement process has a fresh cache and cannot inherit warmth.

`inspect_preparation!` queries that same child without launching a replacement or
extending model lifetime. It checks both retained evidence and every dependent
model. Assigned operations carry the prepared input key; the child rejects a
missing/expired preparation before invoking the operation. The shared
`ScientificResources` agent owner uses this path before work and after successful
completion, and keeps preparation separate from job activity/failure.

Caches bound admitted keys, completed retained bytes and idle lifetime. They
reject new keys when all slots are occupied by builders. Cleanup invalidates
pending identities, so late completion cannot repopulate a cleared cache. The
byte accounting uses `Base.summarysize` for retained values: it is **not** a limit
on transient solver memory, compiled code, or total process memory. Failed-entry
reuse rethrows its stored cause; a newly failed asynchronous builder is surfaced
as Julia's `TaskFailedException`.

Commands and output lines have a finite byte limit. EOF, malformed responses,
deadline/cancellation and failed progress observers cannot leave pending frames
to be mistaken for the next operation's result. Graceful exit calls profile
cleanup; the parent owns hard termination when numerical code cannot cooperate.
Submission itself is also bounded: a full stdin pipe cannot hide cancellation or
the operation deadline. Stop attempts the shutdown command asynchronously, then
TERM/KILL under finite waits and joins owned pipe tasks. Unresolved cleanup retains
the original process handle and throws; callers must not release lease capacity.
This primitive owns one process. Subprocess-tree/container ownership still belongs
to the agent adapter and is not provided by the native process tests.

## Verification and remaining integration

```sh
julia --project=worker/core worker/core/test/runtests.jl
julia --project=worker/core worker/core/test/profile_processes.jl
julia --project=worker worker/test/runtests.jl
```

The profile-process test runs real numerical preparation in separate native Julia
children. It checks cold/warm reuse, process replacement, wrong-profile rejection,
and continued line execution while power flow prepares and after it stops. These
are trusted native tests, not a container isolation or terminal-sandbox claim.

The shared lease-owned scientific manager is implemented and tested with both
these real profiles. Its physical native/container drivers, launch-time policy
enforcement, remote stage/status reporting and targeted broker orchestration are
still being integrated. The public capability endpoint therefore still reports
`assigned_execution=false`, and the runtime controls still report preparation as
unknown. No ready label or runnable agent CLI has been fabricated around these
local process tests.

See [runtime scientific ownership](../../runtime/CONFIGURATION.md#lease-owned-scientific-execution)
for the required driver boundary and additional live-lease numerical test command.
