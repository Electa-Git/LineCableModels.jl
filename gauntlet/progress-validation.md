# Gauntlet progress validation

Validation uses isolated temporary directories and Julia 1.12.7. No existing user
campaign was restarted, resumed, cancelled, or changed. Native boundary tests use
fake GetDP workers and local PSCAD transport/compile fixtures; no real FEM/PSCAD
calculation or additional native warmup/repeat was run.

## Responsibilities and deletions

- `src/progress.jl`, MC and batched traversal report accepted scans with source
  ordering and parent scope identity. Package code has no Gauntlet dependency.
- `gauntlet/progress.jl` synchronously collects outcomes and coarse duration
  evidence, then publishes throttled snapshots. Execution completion, comparison
  verdict and fresh timing eligibility are separate facts. Accepted counts and
  duration observations use their actual completion anchors; duplicate, delayed
  and aborted child observations cannot stretch accepted-scan durations.
- `gauntlet/watch.jl` reads metadata and snapshots in a separate process. It owns
  the six terminal rows, cached spinner and local clocks, input validation, resize,
  color, plain output and cleanup. `--session` pins the invocation across attempts.

The previous calculation-side timer, renderer, logger/output coordination, workload
fingerprint/history learner, FEM mesh-work learner and progress-driven yields were
removed. Required native supervision, native UI cooperation, cancellation, recovery,
scientific outputs and timing records remain owned by execution.

Snapshot schema 2 is disposable. Publication holds one serialization lock through
same-directory `Base.rename`, with a shorter state lock for the payload. Controlled
calls publish suspended state before entering the timer; the post-call transition
records the result and restores observation. Optional IO failure disables that
output; required timing/checkpoint IO still fails normally.

## Regression commands

The default depot is read-only in this environment. Commands used a writable cache:

```bash
export JULIA_DEPOT_PATH=/tmp/lcm-progress-depot:/home/amartins/.julia
julia --project=gauntlet --startup-file=no test/gauntlet/runtests.jl progress_tests
julia --project=gauntlet --startup-file=no test/gauntlet/runtests.jl campaign_tests checkpoint_tests selective_execution_tests campaign_uq_tests
julia --project=gauntlet --startup-file=no test/gauntlet/runtests.jl checkpoint_tests selective_execution_tests
julia --project=test --startup-file=no test/runtests.jl monte_carlo_retries fem_workers pscad/boundaries formulation_overlays feasible_geometry_uq
julia --project=test --startup-file=no test/runtests.jl integration/formulation_grid.jl unit/uq/formulations.jl
python3 -m unittest discover -s test/cli -v
python3 test/gauntlet/progress_terminal.py
```

| Check | Result |
| --- | --- |
| Focused progress and controlled timing | 162 assertions passed |
| Campaign correctness/reuse and UQ | 58 assertions passed |
| Checkpoint and selective execution follow-up | 95 assertions passed |
| Native protocol fixtures, accepted MC retries and geometry UQ | 224 assertions passed |
| Batch equivalence, shared lowering, callbacks and UQ options | 148 assertions passed |
| Python CLI tests | 10 tests passed |
| PTY screen smoke | Passed |

The initial combined campaign run found two test formulations missing the current
report description interface. Their fixture-only `pairs` methods now expose no
native equation children; the checkpoint/selective follow-up passed. A first
boundary command used the wrong project and lacked direct test dependencies; the
reported successful run uses `--project=test`. Visual-tagged tests, including
`formulation_overlays`, remain excluded by the normal runner filter. This report
does not claim the full package, visual suite, or a real native campaign passed.

The PTY smoke launches separate producer and viewer Julia processes, feeds actual
terminal output into a small VT screen emulator, and asserts the resulting screen.
It checks shorter labels without leftovers, a narrow resize and return to six
rows, cached animation during a non-yielding opaque call, truthful controlled-call
suspension, failed comparison with `ETA done`, frozen elapsed/spinner on closure,
and cursor restoration after SIGINT. The CLI enables catchable SIGINT in watch mode
so cleanup runs before exit. Terminal validation used a Linux PTY; Windows
terminal behavior was not exercised.

## Operational overhead

```bash
python3 test/gauntlet/progress_overhead.py /tmp/lcm-progress-overhead.json
```

The reproducible harness runs a fixed-seed owned MC calculation across 48
frequencies, with one Julia compute thread and one BLAS thread. Compilation and the
exact allocation size are warmed before ten interleaved triples: observation off,
publisher only, and publisher with a separately running PTY watcher. Every measured
run has a fresh temporary directory. Watcher startup precedes timing. Numerical
output digests must match across all arms. Publication counts and all raw durations
are retained in the JSON report.

The measured scope is the owned MC operation and its observation/publication
lifecycle. It excludes full campaign reporting/checkpoints and watcher startup; it
does not measure real FEM/PSCAD performance. The approximately one-percent median
overhead target is a measurement target, never a wall-clock regression assertion.
Unrelated machine load remains outside this harness's control.

Measured on 2026-09-13T18:18:08Z, Linux x86_64, Intel Core i9-13950HX,
Julia 1.12.7, one Julia compute thread and one BLAS thread. Calibration
selected **95 accepted scans per run**, each containing 48 frequencies. There
were ten repetitions per arm and 11 snapshot publications in each enabled run
(ordinary updates plus lifecycle boundaries); off produced none. The test processes
were stopped before collecting measurements.

| Arm | Median seconds | IQR seconds | Range seconds | Median paired overhead |
| --- | ---: | ---: | ---: | ---: |
| Off | 6.011370 | 0.003498 | 5.996638–6.027755 | +0.000% |
| Publisher only | 6.007376 | 0.007798 | 5.999755–6.032282 | -0.018% |
| Publisher + external watcher | 6.017132 | 0.014877 | 6.006245–6.046304 | +0.129% |

Paired overhead compares each enabled run with the off run in the same rotated
triple. Its variability was:

| Arm | Paired-overhead IQR | Paired-overhead range |
| --- | ---: | ---: |
| Publisher only | 0.215 percentage points | -0.356% to +0.520% |
| Publisher + external watcher | 0.402 percentage points | -0.169% to +0.619% |

Both measured medians satisfy the approximately 1% target for this warmed owned
workload. The slightly negative publisher-only median is measurement variation,
not evidence that publishing makes computation faster. These results do not imply
zero overhead or predict native-solver/system-load behavior.

Raw durations and paired differences:

| Round | Off (s) | Publisher (s) | Watcher (s) | Publisher Δ | Watcher Δ |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 6.011229447 | 6.011581583 | 6.019257477 | +0.006% | +0.134% |
| 2 | 6.009091646 | 6.006598573 | 6.046303502 | -0.041% | +0.619% |
| 3 | 5.996638280 | 6.006890672 | 6.007748470 | +0.171% | +0.185% |
| 4 | 6.027755402 | 6.006307957 | 6.018238373 | -0.356% | -0.158% |
| 5 | 6.008519174 | 6.004209453 | 6.016026019 | -0.072% | +0.125% |
| 6 | 6.011824573 | 6.015043827 | 6.008905362 | +0.054% | -0.049% |
| 7 | 6.012271549 | 6.007862062 | 6.009053678 | -0.073% | -0.054% |
| 8 | 6.001096420 | 6.032282039 | 6.025340310 | +0.520% | +0.404% |
| 9 | 6.011510644 | 6.025979543 | 6.041872936 | +0.241% | +0.505% |
| 10 | 6.016417145 | 5.999755029 | 6.006245383 | -0.277% | -0.169% |

All 30 numerical output digests matched: `5c6316f5a1c8e2ac1b40640cdbdb7d89323868e74da8839eb7a25e200d3b5311`.

## Use

```bash
./gauntlet/lcm gauntlet run --definition FILE.jl --directory DIR --progress auto
./gauntlet/lcm gauntlet resume --directory DIR --progress auto
./gauntlet/lcm gauntlet status --directory DIR --watch --session SESSION
```

Use the actual quoted watch command printed by run/resume. `plain` adds throttled
execution summaries; `off` removes optional observation while retaining required
timing records. Directory-only `status --directory DIR --watch` remains available
with conservative discovery. Closing the watcher has no effect on execution.
