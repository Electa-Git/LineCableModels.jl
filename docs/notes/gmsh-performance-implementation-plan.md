# Gmsh/GetDP performance implementation plan

Status: implemented and validated. Production batching, isolated frequency
workers, checksummed recovery, cancellation, UI/maps and the four full
101-frequency numerical comparisons are complete. See the
[acceptance report](gmsh-performance-validation.md) and
[machine-readable results](gmsh-performance-validation.json).

The default is two workers with one solver thread each. Full campaign
acceptance used four workers; controlled scaling checks covered one, two,
four and eight workers, with exact serial agreement. This bounds the default
memory requirement while documenting faster settings for the measured cases.

The [diagnosis and experiments](gmsh-performance-refactor.md) provide the
supporting measurements and DSL details. Terminal batching measured about
5.8× acceleration; two concurrent frequency batches measured another 1.98×
relative to sequential batches. These are separate sampled measurements, not
a promised whole-campaign multiplier.

## Target architecture

```text
Julia: resolve physical inputs and prepare/validate all frequency meshes
  -> bounded queue of immutable frequency jobs
     -> independent GetDP processes, each handling one frequency
        -> first requested terminal: assemble and factor
        -> remaining terminals: update RHS and reuse factors
        -> isolated column outputs, optional maps, diagnostics
  -> Julia coordinator: validate, checkpoint, assemble ordered Z/P
  -> existing Engine reductions and condition-checked P-to-Y conversion
```

Each frequency job identifies its mesh, exact frequency/index, material data,
shell radii, requested terminal indices, solver configuration and exclusive
attempt directory. Gmsh and ONELAB remain owned by the coordinator. Workers
receive all numeric inputs explicitly and do not call the shared Gmsh API.

The physical equations, boundary conditions, material resolution, conformal
meshes, terminal ordering and result reductions remain the reference behavior
throughout this work. Existing campaign artifacts remain the comparison
baseline; benchmark runs receive fresh directories and campaign identities.

## Ordered implementation stages

### 1. Define the job and checkpoint contracts

Main locations: `getdp.jl`, run metadata in `compute.jl`, execution options in
`formulations.jl`, and `results.jl`.

- Define frequency jobs and attempt results, including an explicit list of
  missing terminal columns. Separate process launch counts, completed columns
  and completed frequencies.
- Define exclusive attempt paths and independently validated column outputs.
  A column is complete only after Z, P and any requested maps are closed and
  the completion record has been written. Julia validates contents before
  adopting that record.
- Version the execution protocol. Retain existing physical-input, source and
  executable compatibility protections. Keep scheduling metadata such as
  worker count distinct from numerical compatibility in the new protocol.
- Keep a full-assembly, single-column diagnostic route for controlled
  comparisons. Preserve the current retained snapshots as the baseline.

Exit gate: malformed, duplicate, partial and mismatched column artifacts are
rejected; compatible completed columns are identified without mutating them.
Changing scheduling order does not change frequency/terminal identities.

### 2. Integrate terminal batching with one worker

Main locations: maintained `quasi_tem.pro`, `getdp.jl`, and output validation.

- Set runtime source currents for every terminal on every excitation, clearing
  the previous source. Update constraints before RHS generation.
- On the first requested basis, call `Generate` and `Solve`. For subsequent
  bases, call `GenerateRHSGroup[Sys_FEM, Terminals]` and `SolveAgain`.
- Keep one frequency and its mesh fixed for the process lifetime. A resumed
  batch may start at any terminal; first-solve logic must not assume index 1.
- Bind each output and optional map filename/label to its actual basis. Fix
  cross-call append behavior or use basis-specific column operations.
- Establish supported GetDP capabilities explicitly; unsupported binaries
  receive an actionable error or an explicitly recorded diagnostic fallback.

The initial copied-input probes isolated batching from transport. Production
now uses the standalone worker exclusively; the internal full-assembly,
single-column diagnostic controls remain available for calibration.

Exit gate: all requested columns match the unchanged solver at low, middle
and high frequency for detailed/homogenized models and both material profiles.
Confirm one symbolic factorization, one numeric factorization and N solves per
fresh N-terminal frequency batch. Compare reduced results and final Y through
the actual Engine path as well as primitive Z/P.

### 3. Make numerical jobs independent of ONELAB

Main locations: `getdp.jl`, `getdp/model.pro`, and coordinator calls in
`compute.jl`.

- Launch standalone GetDP using Julia `Cmd` arguments, explicit immutable
  input paths and per-child environment settings.
- Give each attempt its own working directory, GetDP `-name` prefix, raw
  outputs, maps and stdout/stderr files.
- Capture exit status and phase timings directly. Publish progress and update
  counters through the coordinator instead of relying on Gmsh logger capture.
- Keep data-only serialization in Julia and solver operations in maintained
  `.pro` files.

Exit gate: standalone and transition-transport runs agree on the same inputs.
An inherited ONELAB value, unrelated Gmsh model, current working directory or
path containing spaces cannot redirect a numerical job's inputs or outputs.

### 4. Add bounded frequency concurrency and recovery

Main locations: `getdp.jl`, execution options, run-state coordination.

- Add configurable worker count and per-solver thread budget. Start validation
  with one and two workers, with one BLAS/OpenMP thread per solver.
- Use a bounded frequency queue. Workers launch/wait for independent GetDP
  processes; one coordinator writes shared metadata and assembles results.
- Adopt complete columns from successful or interrupted attempts after strict
  validation. Retry only missing columns, rebuilding factors in the new
  process. Do not concurrently retry an attempt whose process is still alive.
- Handle worker failures and cancellation explicitly: stop scheduling, reap
  launched processes, retain valid checkpoints and diagnostics, and leave the
  scan incomplete when required results are missing.
- Assemble aggregates once in deterministic frequency/terminal order. Remove
  repeated full-table recounts and cumulative-log rewrites.

Exit gate: one-worker and two-worker results agree; out-of-order completion
does not reorder matrices; injected failures and resumed runs produce the same
complete result without duplicate columns or orphaned processes. Worker count
changes preserve compatible checkpoints under the new protocol.

### 5. Complete UI, maps and public behavior

Main locations: `compute.jl`, `onelab.jl`, result/map validation and
`docs/src/fem.md`.

- Use the same isolated numerical worker for headless and UI execution.
  ONELAB remains the control/progress surface. Gmsh event handling and map
  merging stay with the session owner, outside worker tasks.
- Preserve all nine optional field-map quantities, correct labels and names,
  output basis, tracing, typed errors and session restoration.
- Exercise UI responsiveness and closure/cancellation while jobs are active,
  as well as closure before solving and after successful completion.
- Document batching, worker/thread configuration, process invocation counts,
  recovery and executable requirements.

Exit gate: maps on/off, headless/UI execution, and caller-owned Gmsh sessions
pass their behavioral checks. Existing results and reductions retain their
public meaning.

### 6. Validate the complete campaign and choose defaults

Run relevant Gmsh extension tests, including transport, failure, resume,
enclosure and numerical reference coverage. Add focused tests for the new
batch and concurrent-failure contracts, using small fixtures for repeated
failure injection and real GetDP for numerical/factorization assertions.

Then execute fresh full 101-frequency campaigns for the detailed and
homogenized 18 kV cases and both material profiles, using the retained mesh
snapshots for the numerical comparison. Validate every primitive Z/P slice,
the reduced outputs and final Y against the retained reference. Use tight
roundoff-level tolerances calibrated by unchanged-solver repeats, with an
appropriate absolute tolerance near zero; the existing 10% cross-mesh legacy
threshold is not the acceptance threshold for this refactor.

Benchmark a representative frequency subset with one, two, four and eight
workers as memory permits, keeping solver thread settings explicit. Select a
bounded default from measured throughput and memory, then run the full
campaign with that setting. Record mesh preparation, parsing/preprocessing,
assembly, factorization, back-solves, output and total wall time separately.
Measure one fresh mesh-generation path too, so a cache-only benchmark does not
hide an end-to-end regression. Exercise a small additional multi-terminal
case beyond the 18 kV geometry.

Exit gate: full scans complete with exactly the required columns, maps when
requested, completion markers and checksums; numerical comparisons pass;
profiling verifies one factorization per freshly solved frequency; recovery
tests pass; full-campaign speed and memory use are reported. Broader package
checks required by the touched code must also pass.

## Completion criterion

The work is complete when the production backend performs a full frequency
scan through this worker/coordinator path, preserves the reference numerical
behavior, recovers correctly from interruption, retains the existing UI/map
contracts, and has a reproducible campaign report demonstrating the actual
performance improvement. Probe success alone does not satisfy completion.

Whole-scan GetDP systems, mesh coarsening, alternate physics and more direct
RHS representations remain subsequent investigations. None is required to
deliver the measured terminal-reuse and frequency-parallelism mechanisms.
