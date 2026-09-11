# Campaign progress verification — 2026-09-11

The implementation follows [the execution plan](gauntlet-progress-performance-plan.md).
Observation is task-scoped runtime state, not a numerical declaration field.
No progress hooks were added to `Engine._solve!` or numerical kernels.

## Checks

These are separate runs, with overlapping assertions; their counts are not additive.

| Check | Result |
| --- | --- |
| Focused campaign, UQ, checkpoint, record, local-performance, PSCAD and ownership checks | 701 assertions passed outside the artifact-packaging fixture |
| Artifact packaging and immutable vault protection | 41 passed using a writable temporary Julia depot |
| Final progress, quiet-span, watcher, ETA, terminal, ownership and PSCAD protocol checks | 123 passed |
| Real headless GetDP scan, native timing retention and fresh-process recovery | 77 passed |
| PSCAD worker automation fixture in its own PythonCall environment | 115 passed |
| Updated supervisor parsed on the live Windows station | Passed; no Julia/PSCAD runner found |
| Temporary owned CLI run, followed by terminal-bar recovery | Passed; recovered 3/3 calculation jobs |
| Whitespace and launcher syntax checks | `git diff --check` and `bash -n gauntlet/lcm` passed |

The first packaging check encountered the sandbox's read-only default artifact
depot. Its isolated rerun used `/tmp` for writable artifacts and passed. The
PythonCall worker fixture likewise used a writable temporary depot and the installed
Python 3.11 interpreter; it did not install or launch PSCAD.

The broad sweep was stopped after a profile identified prolonged LLVM optimization
inside the existing 100,000-sample spectral-construction fixture. Relevant campaign,
UQ and adapter checks were then selected explicitly at `-O1`; these are correctness
checks, not performance evidence. The final small checks, CLI checks, real FEM
check and monitoring probe use normal optimization. GUI testing was excluded from
the headless FEM run. The full spectral suite and a new live PSCAD calculation
are not claimed as completed.

## Monitoring-overhead probe

The probe uses `two_insulated_wires`, 2% geometric/layer standard uncertainty,
128 accepted MC trials, seed `0x1234`, and 101 logarithmic frequencies from 0.1 Hz
to 10 MHz. Both paths are warmed before nine alternating-order on/off pairs.
All retained statistics and point seeds are checked for equality.

Monitoring includes in-memory reports, a four-hertz UI timer, plain rendering to
`devnull`, and real throttled snapshot writes in a temporary directory. It does
not measure a particular terminal emulator's output cost. Other validation
processes finish before measurement; this is not a claim of machine isolation.

| Measurement | Monitoring off | Monitoring on |
| --- | ---: | ---: |
| Median compute-call wall time | 15.211512222 s | 15.141887935 s |
| Minimum–maximum | 13.796994529–15.530840358 s | 14.025779666–15.523363497 s |

The ratio of medians gives **−0.4577% observed overhead**. The median target is met
in this probe, but variability is larger than this difference: do not interpret
the negative result as a speedup or a universal sub-1% guarantee. Controlled-sample
silence is also verified structurally, independently of this operational probe.

Raw paired durations in seconds, in execution-pair order:

```text
off: 15.294522851 15.211512222 13.796994529 15.530840358 15.394024941
     15.227891630 14.031984082 15.029100786 15.056628238
on:  15.185158309 15.128365367 14.025779666 14.031703178 15.329634392
     15.523363497 15.199509182 15.141887935 15.125496308
```

Reproduction script for this workspace: `/tmp/lcm_progress_overhead.jl`, run with
`julia --project=gauntlet --compiled-modules=existing --startup-file=no`.
API docstrings explicitly document timing scopes and seconds, following the Julia
docstrings skill; execution, compute-call and native measurements are not conflated.

## Production handoff

The old all-references campaign was interrupted at the user's request. Its saved
campaign directory, FEM runs/meshes and local PSCAD run cache were moved to Trash,
including the small FEM caches created by verification. Definitions, published
reference artifacts, solver installations and the user's REPL were preserved.
The local PSCAD configuration now requests quiet diagnostics independently of
campaign progress. No replacement production campaign was launched.

From the repository root, start the fresh foreground campaign with:

```bash
./gauntlet/lcm gauntlet run \
  --definition gauntlet/benchmark_all_references.jl \
  --configuration gauntlet/local.jl \
  --directory gauntlet/.work/all-references \
  --grid log --count 101 --bounds 0.1,1e7 \
  --on-error continue --progress auto
```

Do not add `--resume`, `--recover-solvers`, backgrounding or log redirection for
this fresh interactive run. The catalogue excludes bare wires from PSCAD.
An empty timing history deliberately shows `estimating`; controlled samples
explicitly pause observation for their complete compute call.
