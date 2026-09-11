# Gauntlet campaign progress and performance: execution plan

Status: implemented on 2026-09-11. See the
[verification and fresh-run handoff](gauntlet-progress-performance-verification.md)
for measured overhead, checks and remaining test-environment limitations.

This plan consolidates the campaign-tracking and performance-accounting decisions.
It does not authorize restarting an existing campaign, changing its scientific
inputs, or adding unrequested expensive solver repetitions.

## 1. Fixed requirements and boundaries

- Default experience: campaign progress visible, routine solver chatter quiet,
  performance recorded with explicit measurement scopes.
- Show completed, failed, running and pending benchmarks; the active case,
  backend, reference/candidate, formulation/parameter point and execution stage;
  Monte Carlo accepted/target trials and retries; external worker/job counts when
  the adapter knows them; elapsed time and qualified remaining-time estimates.
- Do not track or display which frequency a solver is processing. Do not add
  per-frequency hooks to owned numerical kernels.
- Keep campaign progress, diagnostic verbosity and performance measurement
  separate. Progress is not a numerical input or a result-delivery callback.
- During a controlled performance sample, suspend tracker publication, terminal
  redraws and progress-snapshot writes. Show the sample identity before starting
  its timer, and update the display after stopping it. An entire timed Monte Carlo
  call therefore has no live trial updates.
- Preserve scientific results, RNG consumption, trial acceptance, computation
  order, solver batching, recovery behavior and existing user callbacks.
- Do not add fields to problem, formulation or benchmark-declaration structs for
  tracking. No new source-hash gates, scheduler, service, database, event-class
  hierarchy or terminal-dashboard framework.
- Every actual calculation gets scoped operational/source timings. Repeated
  performance measurements run only where declared/requested; do not silently
  enable repeated FEM or PSCAD calculations across the catalogue.

## 2. Current implementation findings

The implementation must address these concrete boundaries rather than treating
all existing elapsed values as equivalent:

- `gauntlet/campaigns.jl::_execute` starts timing before campaign bookkeeping.
  For formulation sweeps its elapsed interval also includes persistence of child
  calculations. The existing `compute_wall` label is consequently too broad.
- `gauntlet/benchmarks.jl::_benchmark_owned` uses `@timed _execute(calculation)`.
  It does not enforce quiet execution. UQ catalogue definitions request these
  repeated checks; FEM and PSCAD catalogue definitions do not.
- `gauntlet/performance.jl::benchmark_local` already provides a warmed
  BenchmarkTools route. The two timing paths need consistent scope and reporting
  rules, not an additional independent timing subsystem.
- GetDP already emits native constraint, assembly, solve and output wall-clock
  durations. FEM workers retain separate process elapsed times and validated
  completion records. Reuse these sources; do not add new numerical probes.
- PSCAD retains remote-machine time around `line.compile()`, explicitly excluding
  output-readiness waiting and transfer, plus broader adapter wall time.
- The Monte Carlo owner already has accepted/attempted/rejected counts. LEP and
  deterministic sweeps already have outer point/formulation boundaries.
- `resume_campaign` currently invokes `run_campaign` for saved definitions in a
  loop. A resumed invocation must have one overall tracker, not repeated 1/1 bars.
- `campaign_status` can inspect saved result identities. A fast watch loop must
  not repeatedly deserialize results or hash their payloads to repaint the screen.

## 3. Minimal architecture

### Runtime reporting in the package

Add `src/progress.jl`, loaded before Engine, ParametricBuilder and UQ.
Use a task-scoped execution-observation context: an optional progress receiver and
a quiet-performance-span flag. Default reporting is a no-op. This context is
runtime-only; never include it in declarations, retained numerical details or
calculation fingerprints.

Use a small number of module-public, non-exported operations following the
repository's ownership/import rules. Do not introduce an abstract event API.
Named tuples are sufficient for reports. Context inheritance must work for nested
calculations and owned tasks, with restoration in `finally`/scoped execution.

Each report identifies the active attempt, benchmark, operand and relevant outer
scope, and carries a stage, unit, absolute completed/total counters and optional
worker/retry counts. A total may be unknown until the owner determines it. Scope
identity includes the relevant point/formulation, so independent calculations do
not overwrite each other.

Counts are absolute, not increment-only notifications. Repeated or coalesced
reports must be harmless. Retain only active scopes and bounded summaries; do not
retain an event or scope for every Monte Carlo trial.

Read the enabled receiver once at an outer execution boundary and pass a concrete
local reporter where necessary. The disabled path must avoid formatting, event
allocation and repeated dynamic lookup inside the trial loop. Numerical kernels
receive no reporter.

### Tracker and display in Gauntlet

Add `gauntlet/progress.jl`, owning the campaign state, ETA calculation, a compact
terminal renderer and the plain-output fallback. One concrete tracker state is
appropriate; no generic graph or subscription framework is needed.

Gauntlet owns benchmark/operand/stage lifecycles. Domain owners report their own
work; they do not print progress bars or import Gauntlet. Existing `on_result`
callbacks retain their current purpose and behavior.

Use a single synchronized state update path for owned concurrent tasks. Updates
replace current counters and do not append to an unbounded channel. Only the
tracker renders or writes its optional snapshot. Snapshot data is copied under
the state lock; do not hold that lock during terminal or filesystem IO.

## 4. Counting and lifecycle semantics

- The run denominator is the selected benchmark definitions for that invocation,
  including all selected definitions on resume. The watch command also exposes
  the persisted campaign inventory. Label a subset run as selected work.
- A benchmark is not finished until required calculations, validation, reports,
  persistence and any declared performance checks have finished.
- Track successful, failed, skipped and interrupted outcomes separately. A final
  100% finished bar is not a claim that every benchmark passed.
- Keep calculation jobs, external worker jobs and MC trials as separately labelled
  units. Never add them together into one denominator.
- Derive outer job counts from existing declarations and traversal cardinalities,
  without materializing random draws or eagerly expanding an enormous MC queue.
  Use unknown totals where the execution owner has not resolved the work yet.
- Successful reuse is a subset of completed work. Count a result as reused only
  after existing recovery validation accepts it. Preserve original timing and
  provenance; do not call a cache hit a zero-second solver measurement.
- Distinguish whole-result recovery, partial native recovery and intentional
  within-solve factorization sharing. The latter is part of the numerical method,
  not a reason to discard its performance measurement.
- Failed MC draws increment attempts/rejections, not accepted trials. An inferred
  MC target becomes known after the existing target-selection logic runs.
- Failed dependencies make downstream work skipped, not permanently pending.
  Interruptions close the current tracking scope without claiming success or
  changing the package's existing cancellation/recovery policy.
- A stale heartbeat means no recent observation. It does not by itself prove
  solver failure or authorize cancellation.

## 5. Timing contract and controlled samples

Keep these measurements distinct in retained records and reports:

| Scope | Meaning |
| --- | --- |
| Campaign/session wall | End-to-end invocation duration, including preparation, IO, reporting and monitoring. Previous downtime is not new execution time. |
| Execution/orchestration wall | Actual workflow elapsed around one operand, with its bookkeeping and persistence stages identified. |
| Compute-call wall | Time around the actual `compute` call, excluding Gauntlet bookkeeping and persistence. Backend-internal preparation, IO and any callbacks executed by that call remain included; record the callback/diagnostic policy. |
| Backend/source timings | Native GetDP phase durations, remote PSCAD compile duration and worker-process elapsed, each retaining its source and exclusions. |
| Controlled compute samples | Repeated, appropriately warmed compute-call timings with diagnostic/progress output suppressed, labelled with sample count, execution settings and environment. |

Place monotonic timers at these boundaries and publish timing records after the
measured work. For sequential formulation sweeps, aggregate child compute spans
without including their persistence or double-counting parent elapsed time.
Parallel worker durations are accumulated worker time, never summed and labelled
as elapsed solve time. Keep output, startup, transfer and preparation distinct
from native numerical solve phases where the source provides that distinction.

The observer is a consumer of timing records, never their authoritative source.
ETA estimates are not performance measurements and must not enter speedup tests.
Julia allocation measurements describe Julia allocations, not total native-worker
memory usage.

Controlled performance execution must:

1. Resolve declared benchmark inputs and prepare the timed call outside the
   measurement boundary. Sampling, realization/reconstruction, retries and
   aggregation belonging to a Monte Carlo computation stay inside its measured
   call; do not time only its accepted inner solves and call that MC performance.
2. Warm the relevant owned Julia path outside measurement when required. Record
   the warmup policy; do not launch native warmup/repeat runs implicitly.
3. Disable the progress receiver and pause the renderer/snapshot producer for the
   whole measured call. Synchronize with any in-flight renderer before timing.
4. Honor task-scoped diagnostic silence in the owned loggers and external solver
   launchers; an outer `NullLogger` alone is insufficient because compute methods
   currently install their own console loggers. Do not use process-global logging
   switches or global stdout redirection to silence unrelated work.
5. Keep progress and ordinary user `on_result` delivery outside controlled samples.
   Preserve those callbacks in the normal correctness/execution pass. Suppression
   in the separate performance pass is explicit and recorded, not silently counted
   as timing the callback-bearing workload.
6. Retain required solver output/checkpoint IO as part of the measured algorithm;
   silence optional chatter, not numerical outputs. Keep settings and timing
   scope explicit so quiet and verbose runs are not silently conflated.
7. Record samples, then refresh the tracker outside the timed region. Errors and
   interrupts still propagate; a quiet span must not swallow calculation failures.

Use existing `ParametricProblem.options` ownership when preparing UQ execution;
do not add unsupported keyword forwarding to UQ `compute` methods or mutate saved
declarations. Any runtime diagnostic override must be explicit in the measurement
record. Existing correctness/reporting inputs remain unchanged.

The comparability decision must consider scope, execution settings, reuse,
environment and measurement policy, not only the current coverage/allocation-log
check. Report limited evidence honestly when only one timing sample is available.
No claim of zero overhead or isolation from unrelated system load is permitted.

## 6. Instrumentation by owner

### Gauntlet and owned deterministic calculations

Instrument campaign, operand and existing outer point/formulation boundaries, plus
preparation, validation, report, save and performance-sample stages. Do not add
hooks in `Engine._solve!` or matrix/earth-return kernels. Scalar opaque calculations
show their stage and elapsed time until returning.

### Monte Carlo and LEP

Outside controlled samples, expose the MC loop's existing counts after accepted
trials or rejected attempts, using a cheap in-memory update. Rendering and disk IO
are not performed by that loop. Preserve the exact RNG and acceptance sequence.
LEP reports outer traversal/computation stages, not a fabricated trial counter.
Aggregation has its own stage so finishing the trial count does not imply that
the result is already saved.

### FEM/GetDP

Use `_start_worker!`, `_finish_worker!`, `_record_progress!` and existing run-state
transitions. Show active, completed, queued and recovered work using the adapter's
actual worker/validation state. Do not expose frequency identifiers. Consume
existing native timing files and checkpoint records; do not re-read or hash every
artifact on each UI refresh. Track meshing, solve, parse/validation and finalization
as distinct stages without changing worker scheduling or factorization reuse.

### PSCAD

Report staging, launch, loading/configuration, compile, output-readiness,
transfer and validation stages where the owner can observe them. Use the existing
remote supervisor's process checks for liveness. No frequency-progress investigation
or scan splitting is needed for this scope.

Carry small, explicitly identified machine-readable progress records over the
existing remote transport, independently of human log verbosity. Do not parse
free-form solver messages. The supervisor forwards observations even when ordinary
console streaming is disabled. A runner heartbeat is not a statement that the
native solver has advanced. During controlled samples, suppress monitoring
publication while preserving timeout, cancellation and required output handling.

Preserve `timing.txt` and its current compile-only meaning; use a remote monotonic
duration measurement. Do not substitute client transport time for source time.

## 7. ETA policy

- Maintain separate estimates for the active calculation and the remaining run.
- Base MC estimates on recent accepted-trial throughput, including rejection cost;
  separate setup from the steady-state estimate.
- Base FEM estimates on observed completed-job throughput at the actual worker
  count and comparable mesh/workload settings, with startup/meshing kept separate.
  Known job-cost differences and partial recovery must not be treated as equal-cost
  completion percentages.
- Use comparable previous PSCAD whole-run/stage observations when its active stage
  is opaque. Lack of a new observation must not count as progress.
- Seed estimates from existing compatible timing records, read once and indexed
  in memory. Group by backend, workload/case, relevant execution settings and
  environment. These are advisory similarities, not new numerical reuse gates.
- For the current sequential campaign scheduler, combine remaining active work,
  pending calculations and required finalization/performance work. Do not sum
  concurrent worker durations as though those workers execute sequentially.
- Do not infer FEM/UQ cost from the PSCAD-only prefix of the all-reference queue.
  Show `estimating` and identify uncovered work when history is insufficient.
- Show approximate rounded estimates or an empirically supported range, never
  unsupported statistical confidence. Estimates may rise when execution slows;
  do not clamp an overdue running stage to zero remaining time.
- Recovered work is removed from remaining solve work only after recovery is
  confirmed. Recovery validation and finalization still take wall time.

## 8. UI, watch mode and compatibility

Implemented interfaces:

```julia
run_campaign(directory, definitions; progress=:auto, ...)
resume_campaign(directory; progress=:auto, ...)
```

```bash
./gauntlet/lcm gauntlet run ... --progress auto
./gauntlet/lcm gauntlet resume ... --progress auto
./gauntlet/lcm gauntlet status --directory DIR --watch
```

Support `auto`, `plain` and `off`. Auto uses a compact updating display on a
supported terminal and plain output otherwise. Use stderr for progress. Respect
terminal width and color support; retain a final readable summary.

Illustrative ordinary-run display:

```text
Benchmarks [████░░░░░░░░░░░░░░░░] 12/59 finished
Successful 11 · Failed 1 · Running 1 · Pending 46
Active 13/59 · NA2XS2Y trefoil · reference / Monte Carlo
Trials 96/512 accepted · Attempts 99 · Rejected 3
Elapsed 00:42:18 · Active ETA ~00:08:40 · Campaign ETA: estimating
```

During a controlled sample, replace the active detail with `Performance sample
2/3 — display paused until sample finishes`. Do not keep drawing a fake live bar.

Maximum ordinary terminal refresh rate: 4 Hz. Plain output: stage/outcome changes
and a throttled periodic update, approximately every 10 seconds. Coalesce rapid
stage transitions. Warnings/errors must remain readable; routine native chatter
should not fight the display. Preserve explicit user diagnostic choices and label
their effect on operational timings.

Write at most one small progress snapshot per second and at important transitions,
under the existing campaign sessions directory, keyed by execution session. Include
active attempt identities and the timestamp; keep the payload bounded. This is
disposable UI state, not another checkpoint or source of numerical truth.

`status --watch` reads the lightweight snapshot and existing campaign state, with
ownership/staleness checks. It must not deserialize executable declarations,
recalculate reports, hash numerical results repeatedly, or launch solvers. Legacy
runs without snapshots still get coarse campaign status. Frozen snapshots during
timing spans are explicitly marked as such. Multiple invocation snapshots cannot
overwrite one another or masquerade as a newer attempt.

`progress=off` disables progress reporting/rendering/snapshot production, not
performance records or authoritative campaign state. Optional display/snapshot IO
failures disable that output with one warning; required result persistence failures
still fail normally. Always restore the terminal and scoped context on exit.

Add new timing metadata compatibly. Preserve legacy fields/readability but label
their actual legacy elapsed scope; do not reinterpret old wall times as new clean
compute samples or rewrite old results. Ordinary display switches do not change
scientific identities, native-recovery compatibility or saved declaration layouts.

## 9. Implementation order and file ownership

1. **Timing contract and regression fixtures.** Define additive timing records and
   tests for scopes, bookkeeping exclusion, reuse and legacy reads. Main files:
   `gauntlet/campaigns.jl`, `benchmarks.jl`, `performance.jl`, `records.jl`,
   `comparisons/saved.jl` and existing campaign/retention tests.
2. **Runtime reporting and tracker lifecycle.** Add the small shared package file
   and `gauntlet/progress.jl`; wire run/resume/operand/traversal scopes. Keep outer
   callbacks and declared calculation inputs unchanged.
3. **Quiet performance spans.** Tighten measured boundaries, warmup policy and
   scoped quiet handling in owned loggers/UQ/backend launchers. Cover complete
   sample isolation before any tracker is presented as performance-safe.
4. **Backend/MC reporting.** Connect existing MC counters, FEM worker transitions
   and PSCAD supervisor stages/liveness. Normalize retained native timing summaries
   without adding frequency or numerical-kernel instrumentation.
5. **Display, snapshots and CLI.** Implement terminal/plain/off behavior, unified
   resume totals, watch mode, safe output lifecycle and failure handling.
6. **ETA and performance presentation.** Add bounded throughput histories,
   source-record estimates, explicit unknowns and scope-correct report summaries.
7. **Verification and documentation.** Run the acceptance checks below and update
   `gauntlet/README.md`, the Gauntlet user documentation and relevant Julia
   docstrings. Provide tested commands and limitations before campaign handoff.

Do not modify numerical formulas, physical benchmark options, frequency grids,
comparison bands, MC seeds/trial targets, solver concurrency or native recovery
algorithms as part of this work. Default catalogue console settings may be made
quiet; explicitly supplied diagnostic settings remain respected outside controlled
samples.

## 10. Acceptance tests and handoff

### Correctness and ownership

- Progress on/off gives identical numerical results, seeds, trial counts, rejection
  records and scientific inputs. Existing `on_result` callbacks fire with their
  existing arguments in normal execution, without duplication or replacement.
- Cover deterministic scalar/sweeps, LEP, MC with retries and inferred trial count,
  unknown-duration backends, native recovery, failed dependencies and Ctrl-C.
- Resuming several benchmarks has one stable run denominator. Restoring a complete
  result performs no extra solve; partial FEM recovery does not double-count jobs.
- Exercise duplicate/coalesced observations and overlapping task scopes. Memory
  remains bounded with many trials and long campaigns.
- Current ownership, ExplicitImports and semantic-economy tests remain green.
- Confirm no frequency-progress call has been added to owned numerical loops.

### Timing and overhead

- Use deterministic fake timing/IO fixtures to prove that bookkeeping, report
  generation, callbacks and tracker output are outside controlled samples.
- Test quiet spans against actual owned logger installation, nested UQ compute,
  native launcher settings, renderer activity and snapshot writes; restore state
  after errors and interrupts. Keep explicit measurement-policy metadata.
- Check native timing-file values against retained summaries; verify units/scopes,
  parallel worker accounting and distinction between reuse and new execution.
- Compare warmed progress-on/off operational runs using repeated measurements;
  record overhead and variability. Target no more than 1% median overhead on a
  representative long owned/MC workload, with enough repetitions to distinguish
  the effect from noise. If exceeded, reduce reporting work before accepting.
- Verify controlled sample execution has no reporter invocations, event allocation,
  terminal redraws or optional progress-snapshot writes. This is a structural test
  in addition to timing measurements, not a claim of a noise-free machine.
- Test ETA with a fake monotonic clock, heterogeneous backend costs, recovery,
  worker concurrency, rejected trials, missing history and overdue work. No
  negative, infinite or NaN displayed ETA; no invented whole-campaign estimate.

### UI, persistence and integration

- Test terminal and redirected output, narrow terminals, warning interleaving,
  broken pipes, failed optional snapshot writes and terminal cleanup.
- Test watch mode against current, legacy, stale, interrupted and frozen timing
  snapshots. Confirm it does not deserialize models or invoke numerical work.
- Verify old campaign/results remain readable and progress switches do not alter
  declaration layouts or numerical reuse decisions.
- Run focused current-tree Gauntlet, UQ, logging, persistence, Gmsh and PSCAD adapter
  suites. Use temporary campaign directories and existing solver fixtures first.
- Once implementation is ready and authorized, verify a small owned/MC/LEP run,
  one representative real FEM run and one PSCAD run in a separate smoke directory.
  Never include bare wires in PSCAD; do not disturb an active user campaign.

Handoff includes the test results, measured monitoring overhead, examples for
CLI/REPL/watch usage and remaining limitations. Do not claim the full all-reference
campaign ran unless it actually did. Existing campaign restart/resume is a separate
execution decision, made only after inspecting its current ownership and saved
declaration compatibility.
