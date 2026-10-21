# Test suite

Tests check the current implementation and architecture, with a **95% production
line-coverage gate**. Scientific acceptance is a research judgment. Float32 inputs
must work without type-induced crashes. There is no extra Float32 accuracy promise.
See the [testing requirements](../docs/src/developers.md#testing-requirements).

## While coding

Run from the repository root with Julia 1.12. Resolve and instantiate the environment before first
use and after dependency changes:

```sh
julia --project=test -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
```

Use the corresponding project for additional environments. `resolve` refreshes
local manifests after changes to the developed package. `instantiate` installs
and precompiles that resolved graph. Finish those operations before starting tests.
Ordinary tests share the root workspace `Manifest.toml`. There is no separate
`test/Manifest.toml`. Do not run these environment operations concurrently or use `--compiled-modules=no` to
work around stale manifests. That flag is retained only for the existing native
and visual commands below.

```sh
julia --project=test test/runtests.jl
julia --project=test test/runtests.jl unit/importexport/atp
julia --project=test test/runtests.jl integration/feasible_geometry_uq
julia --project=test test/runtests.jl tag:quality
julia --project=test test/runtests.jl changed:HEAD
```

Use a file selection while editing. The unfiltered command is the full sweep.
The unfiltered Julia command runs ordinary unit, integration and non-native
extension checks, including small UQ sampling, aggregation and linear propagation.
`tag:quality` checks architecture, formula ownership, descriptions and the actual
external interfaces of loaded adapters.

Selectors match substrings in filenames or item names, or exact `tag:<name>`. Multiple names are
alternatives. Multiple tags are alternatives. Both groups must match when combined.
Append `--list` to inspect selection without executing bodies. Empty selections
and unknown options fail. The runner prints item starts, selected files and items,
elapsed time and completion or failure. A started item is not necessarily completed.

`--workers=N` runs the selection in N new Julia processes. The default is a quarter of
the CPU threads, at most 8. Under CI, which sets `CI`, the default is 1. With one worker,
the selection runs in the current process. `--list` always lists in the current process.
The run fails when any process fails.

With more workers, each worker takes whole test files. The runner packs the files by
the item durations that earlier runs measured, kept in `.linecablemodels/test/durations.tsv`,
or by item count without them. The `quality` items run in a separate new process that
starts with the workers, or in two when one would take longer than the longest worker.
Each process has its own `TMPDIR` and its own log under
`.linecablemodels/test/<run>/`. The console shows each item start with its process, then a
summary line per process and the total.

`changed:REF` selects the items that the changes since the git revision `REF` can
affect, untracked files included. It takes no other selector. The runner uses the
test taxonomy below to choose them:

- every item whose owner tag is the owner of a changed file under `src/` or `ext/`,
  or a later owner. `slow` items are left out.
- every item of a changed test file, `slow` or not, and every item that uses a setup
  defined in a changed test file.
- after a change under `test/support/`, every item that is not `slow`, and every item
  under `test/unit/core/`, which checks the runner itself.
- every `quality` item.

Items with an environment tag stay out, as in an ordinary run. With one worker, the
selected items outside `quality` run in the current process. The `quality` items then run
in a new process, exactly as `tag:quality` runs alone, because the architecture guards
read the live method table that earlier items can extend. With more workers, that
process starts beside the workers. The run fails when any part fails. `--list` lists
both parts.

## Test taxonomy

Each item outside `quality` has exactly one owner tag. The owner tag is the latest layer,
in load order, whose code the item exercises. A change in one layer can affect only that
layer and the layers loaded after it. `changed:REF` relies on this order. The owner tags
in load order are `units`, `commons`, `materials`, `earth`, `datamodel`, `engine`,
`parametric`, `uq`, `report`, `importexport` and `pscad`. The extension owners
follow them: `measurements`, `distributions`, `xlsx` and `makie`. `quality`
items have no owner tag.

`test/support/taxonomy.jl` defines the owner tags, the package module that each owner
covers and the kind tags. The runner, the architecture guards and the tools read it.
It also gives each source file the owner of its load position. A file that
`src/LineCableModels.jl` includes directly takes the owner of the last module included
before it, or the first owner, so `src/docstrings.jl` belongs to `units`. TextDisplay
and PlotBuilder belong to `commons`. The files of a child module take the owner of its
parent: the formula families and ModalAnalysis belong to `engine`.

Each item also has at least one kind tag: `unit`, `integration`, `extension`,
`visual`, `quality` or `aqua`. The environment tags `visual`, `pscad_native`,
`core_only` and `aqua` keep items out of ordinary runs. `slow` marks
items that took more than 10 s in the last serial full run (`--workers=1`).

The quality items check the taxonomy. An item outside `quality` has one owner tag and a
`quality` item has none. Each item has a kind tag and only known tags. An item under
`test/unit/<dir>/` has the owner of `src/<dir>` or a later one.

Give a new item the owner of the latest layer it reaches. To check the choice, record
the code that each item runs and compare it with the tags:

```sh
julia --project=test --code-coverage=@$PWD --code-coverage=/tmp/trace.info test/tools/owners.jl record DIR
julia --project=test test/tools/owners.jl report DIR
```

The report lists every item whose owner tag is earlier than the code it runs.

To review the `slow` tags, run the full ordinary suite with `--workers=1` and give its log,
and the log of each environment run, to `test/tools/durations.jl`. Worker logs cannot
decide the tags. Each fresh worker compiles again what one process compiles once. An item
takes 1.5 to 3 times longer there, and the factor changes with the items before it. The
tool derives the duration of each item from the `Starting` lines.

It lists the items above 10 s without `slow`, the `slow` items below 5 s and the minutes
per owner tag. It reports these findings and leaves the tags unchanged.

```sh
julia --project=test test/tools/durations.jl .linecablemodels/test/RUN/worker-*.log visual.log
```

## Cadence

During a change, run `changed:REF` from the commit that the step starts from. At the end
of a track, before a push, run the full suite including `slow` items. Also run the
extension and visual environments that the change touches, the documentation build, and
Vale and cspell with the versions that CI uses. Finally, run `test/tools/durations.jl` on
the new local logs and the equivalence check that the
[developer guide](../docs/src/developers.md#preservation-locks) describes. The timing
comparison runs in CI and is advisory, with a warning as its strongest outcome. The end of
a track runs every deferred test.

## Advisory source diagnostics

`tag:quality` also runs Fatou and three repository-owned ReLint rules against every
Julia source file under `src/` and `ext/`, including PSCAD. These items read source
without importing LineCableModels or starting external solvers. Ordinary runs
exclude them, and `--list` executes neither scanner.

Resolve and instantiate the test project to install ReLint at revision
`61079427d9fc91cb3ead8f5b287aa3ae9b2269cb` (0.9.0) and Argus 0.4.2. Install
Fatou **0.22.0** separately and place it on `PATH`. For Linux or macOS, use the
same pinned installer as CI:

```sh
curl --fail --location \
  https://raw.githubusercontent.com/jolars/fatou-action/6b6a66a1f28fce0d1ff31048c609e63b7006136e/scripts/install-fatou.sh \
  --output /tmp/install-lcm-fatou.sh
FATOU_VERSION=v0.22.0 FATOU_VERIFY_CHECKSUM=true \
  FATOU_INSTALL_DIR="$HOME/.local/bin" sh /tmp/install-lcm-fatou.sh
export PATH="$HOME/.local/bin:$PATH"
fatou --version

julia --project=test test/runtests.jl tag:quality --list
julia --project=test test/runtests.jl quality/fatou.jl quality/relint.jl
julia --project=test test/runtests.jl tag:quality
```

Windows binaries are available from the
[Fatou 0.22.0 release](https://github.com/jolars/fatou/releases/tag/v0.22.0).
Tests never install tools. A wrong Fatou version fails the item. A missing Fatou
fails the item in CI, which sets `CI=true`. A local run without Fatou records the item
as skipped and prints a message.

All selected source findings are **advisory**, regardless of their count.
`fatou.toml` selects 15 rules explicitly and sets their severity to `warning`.
ReLint runs only `runtime-eval`, `constant-catch-result`, and `private-forwarder`.
The scanner can report typed forwarding methods as candidates without determining
whether their dispatch purpose is useful. Inert quotations and the documented
syntactic exclusions are exercised by controls in each test item.

Invalid configuration, unreadable source, parse failures, malformed output,
unexpected termination, and failed controls fail the quality job. Fatou's
findings-only exit status is accepted after validating its diagnostics. A passing
item establishes that its controls and scan completed. It does not establish the
absence of findings or certify architecture or scientific behavior.

Each scanner prints every finding and per-rule counts. Fatou also prints its
native JSON, including parse diagnostics. CI retains the complete quality log in
the `quality-test-report` artifact, including failed runs that produced a log.
Source-file counts do not contribute to the production line-coverage requirement.
No automatic repairs, formatting, suppression additions, or promotion to blocking
source rules are part of these checks.

## Additional execution environments

| Command | What it executes |
| --- | --- |
| `julia --project=test/core test/runtests.jl tag:core_only` | Optional-extension behavior with core dependencies only. |
| `LINECABLEMODELS_TEST_PLOTTING=true julia --project=test/visual --compiled-modules=no test/runtests.jl tag:visual` | Existing rendering, plot lifecycle, layout and output checks. No new calibration. |
| `julia --project=test test/runtests.jl tag:pscad_native` | Real PSCAD acceptance. Set `LINECABLEMODELS_PSCAD_CONFIG` to your private TOML filename. Requires the configured live station. |
| `julia --project=test test/runtests.jl tag:aqua` | Aqua in a fresh Julia process. |
| `julia --project=docs docs/doctest.jl` and `julia --project=docs docs/make.jl` | Docstrings and documentation build. Instantiate with `docs/instantiate.jl`. |
| `julia --project=src/pscad/remote test/integration/pscad/remote_runner_protocol.jl` | Remote-runner protocol against a mocked Python automation interface. Does not execute PSCAD or validate solver numerics. |

Ordinary exclusions are `quality`, `aqua`, `visual`, `core_only` and `pscad_native`.
These tags describe purpose and environment, never outcome.
PSCAD native tests require the explicit `tag:pscad_native` selector. A filename or item-name
search alone cannot launch a station. Ordinary local PSCAD tests use Julia protocol
fixtures. The dedicated remote-runner protocol test mocks only the Python automation
surface needed to verify runner orchestration. It does not imitate PSCAD numerics.
Native coverage gaps remain visible in the unchanged production coverage inventory
and 95% gate.
There is no required scientific-study job or `scientific` selection.

## Coverage and release verification

The [existing CI workflow](../.github/workflows/CI.yml) runs ordinary tests on
Julia 1.12 and prereleases. It also checks the environments above, Cairo/GL/WGL
activation, clean installation, documentation and merged coverage. Release
verification adds execution environments. It does not certify model physics.

For coverage, use one Julia/BLAS thread and remove old traces first:

```sh
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
julia --project=test/coverage test/coverage.jl clean
julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'
```

Run the additional environments above with `--code-coverage=@.`. For activation,
use `julia --project=test/visual --code-coverage=@. test/runtests.jl "loaded extension activation"`
for Cairo. CI contains the temporary-environment commands for GLMakie and WGLMakie.
Then merge and enforce the gate:

```sh
julia --project=test/coverage test/coverage.jl check
```

This inventories every production Julia file under `src/` and `ext/`, writes
`lcov.info` and fails below 95%. The check removes traces afterward, including on failure.
Retain the report. Missing executions and failing tests remain visible even if
the line percentage passes.

## Measured execution and remaining work

The [execution report](../local/validation-refoundation/2026-09-15/execution/report.md)
records scopes and measured wall times on the shared local machine.
The ATP file selection took 41 s. Core-only environment setup plus six checks took
74 s. Quality took 279 s. The visual selection took about 25 minutes. These include compilation and are not isolated CI timing
predictions.

The archived runs describe an earlier solver inventory and are not validation
of the current QuadGK-only spectral path. The Julia 1.12 full run reached its
90-minute limit. Final coverage is **19,609/20,582 production lines (95.27%)**, passing the unchanged
95% gate. Three earth-return items remain unfinished on Julia 1.12. They completed
on 1.13.0-rc4. The report identifies these exact items. Coverage does not close them.

The earlier FEM/analytical mutual-admittance difference, transformed-exterior
quadrature sensitivity and unresolved derivative-reference study remain in the
[dated evidence](../local/validation-refoundation/2026-09-15/test-system-audit.md).
These research observations require scientific review. They do not establish a code failure or approve a release.
