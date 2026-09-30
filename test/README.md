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

## Additional execution environments

| Command | What it executes |
| --- | --- |
| `DISPLAY= julia --project=test --compiled-modules=no test/runtests.jl tag:fem_numerical` | Native GetDP/Gmsh execution, extraction, material transport, reductions, failure and resume. The UI item is skipped here. |
| `julia --project=test test/runtests.jl extensions/fem_ui.jl` | Four UI lifecycle scenarios in fresh processes. Requires an accessible `DISPLAY`. CI uses `xvfb-run -a`. |
| `julia --project=test/core test/runtests.jl tag:core_only` | Optional-extension behavior with core dependencies only. |
| `LINECABLEMODELS_TEST_PLOTTING=true julia --project=test/visual --compiled-modules=no test/runtests.jl tag:visual` | Existing rendering, plot lifecycle, layout and output checks. No new calibration. |
| `julia --project=test test/runtests.jl tag:pscad_native` | Real PSCAD acceptance. Set `LINECABLEMODELS_PSCAD_CONFIG` to your private TOML filename. Requires the configured live station. |
| `julia --project=test test/runtests.jl tag:aqua` | Aqua in a fresh Julia process. |
| `julia --project=docs docs/doctest.jl` and `julia --project=docs docs/make.jl` | Docstrings and documentation build. Instantiate with `docs/instantiate.jl`. |
| `julia --project=ext/LineCableModelsPSCADExt/remote test/integration/pscad/remote_runner_protocol.jl` | Remote-runner protocol against a mocked Python automation interface. Does not execute PSCAD or validate solver numerics. |

Ordinary exclusions are `quality`, `aqua`, `visual`, `core_only`, `fem_numerical` and `pscad_native`.
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
Complete the display-dependent UI run as well. Then merge and enforce the gate:

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
74 s. Quality took 279 s. Native FEM and the existing visual selection took
about 25 minutes each. These include compilation and are not isolated CI timing
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
