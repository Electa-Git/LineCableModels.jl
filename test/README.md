# Test suite

Tests check the current implementation and architecture, with a **95% production
line-coverage gate**. Scientific acceptance is a research judgment. Float32 inputs
must work without type-induced crashes; there is no extra Float32 accuracy promise.
See the [testing policy](../docs/src/developers.md#testing-policy).

## While coding

Run from the repository root with Julia 1.12. Prepare the environment before first
use and after dependency changes:

```sh
julia --project=test -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
```

Use the corresponding project for additional environments. `resolve` refreshes
local manifests after changes to the developed package; `instantiate` installs
and precompiles that resolved graph. Finish preparation before starting tests.
Ordinary tests share the root workspace `Manifest.toml`; there is no separate
`test/Manifest.toml`. Do not launch overlapping cold preparations or use `--compiled-modules=no` to
work around stale manifests. That flag is retained only for the existing native
and visual commands below.

```sh
julia --project=test test/runtests.jl
julia --project=test test/runtests.jl unit/importexport/atp
julia --project=test test/runtests.jl integration/feasible_geometry_uq
julia --project=test test/runtests.jl tag:quality
```

Use a file selection while editing; the unfiltered command is the full sweep.
The unfiltered Julia command runs ordinary unit, integration and non-native
extension checks, including small UQ sampling/aggregation and linear propagation.
`tag:quality` checks architecture, formula ownership, descriptions and the actual
external interfaces of loaded adapters.

Selectors match file/name substrings or exact `tag:<name>`. Multiple names are
alternatives; multiple tags are alternatives; both groups must match when combined.
Append `--list` to inspect selection without executing bodies. Empty selections
and unknown options fail. The runner prints item starts, selected files/items,
elapsed time and completion/failure. A started item is not necessarily completed.

## Additional execution environments

| Command | What it executes |
| --- | --- |
| `DISPLAY= julia --project=test --compiled-modules=no test/runtests.jl tag:fem_numerical` | Native GetDP/Gmsh execution, extraction, material transport, reductions, failure and resume. |
| `LINECABLEMODELS_TEST_PLOTTING=true julia --project=test test/runtests.jl integration/native_makie_fem.jl` | Saved mesh and field inspection through Makie without GetDP or a display. |
| `julia --project=test/core test/runtests.jl tag:core_only` | Optional-extension boundaries with core dependencies only. |
| `LINECABLEMODELS_TEST_PLOTTING=true julia --project=test/visual --compiled-modules=no test/runtests.jl tag:visual` | Existing rendering, plot lifecycle, layout and output contracts; no new calibration. |
| `julia --project=test test/runtests.jl tag:pscad_native` | Real PSCAD acceptance. Set `LINECABLEMODELS_PSCAD_CONFIG` to your private TOML filename; requires the configured live station. |
| `julia --project=test test/runtests.jl tag:aqua` | Aqua in a fresh Julia process. |
| `julia --project=docs docs/doctest.jl` and `julia --project=docs docs/make.jl` | Docstrings and documentation build; instantiate with `docs/instantiate.jl`. |
| `julia --project=ext/LineCableModelsPSCADExt/remote test/integration/pscad/remote_runner_boundary.jl` | Remote-runner protocol against a mocked Python automation boundary; does not execute PSCAD or validate solver numerics. |

Ordinary exclusions are `quality`, `aqua`, `visual`, `core_only`, `fem_numerical` and `pscad_native`.

The detached ONELAB numerical tests invoke native Gmsh and GetDP executables.
They require no interpreter or extra package environment in the detached bundle.
Julia is used only by the repository test harness to export and check results.

## Coverage and release verification

The [existing CI workflow](../.github/workflows/CI.yml) owns the complete recipe:
Julia 1.12 and prerelease ordinary runs, all environments above, Cairo/GL/WGL
activation, clean installation, documentation and merged coverage. Release
verification adds execution environments; it does not certify model physics.

For coverage, use one Julia/BLAS thread and remove old traces first:

```sh
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
julia --project=test/coverage test/coverage.jl clean
julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'
```

Run the additional environments above with `--code-coverage=@.`. For activation,
use `julia --project=test/visual --code-coverage=@. test/runtests.jl "loaded extension activation"`
for Cairo; CI contains the temporary-environment commands for GLMakie and WGLMakie.
Complete the independent FEM inspection run as well. Then merge and enforce the gate:

```sh
julia --project=test/coverage test/coverage.jl check
```

This inventories every production Julia file under `src/` and `ext/`, writes
`lcov.info`, and fails below 95%. The check removes traces afterward, including on failure;
retain the report. Missing executions and failing tests remain visible even if
the line percentage passes.

## Measured execution and remaining work

The [execution report](../local/validation-refoundation/2026-09-15/execution/report.md)
records scopes and measured wall times on the shared local machine.
The ATP file selection took 41 s. Core-only preparation plus six checks took
74 s; quality took 279 s; native FEM and the existing visual selection took
about 25 minutes each. These include compilation and are not isolated CI timing
predictions.

The archived runs describe an earlier solver inventory and are not validation
of the current QuadGK-only spectral path. The Julia 1.12 full run reached its
90-minute limit. Final coverage is **19,609/20,582 production lines (95.27%)**, passing the unchanged
95% gate. Three earth-return items remain unfinished on Julia 1.12; they completed
on 1.13.0-rc4. The report identifies these exact items. Coverage does not close them.

The earlier FEM/analytical mutual-admittance difference, transformed-exterior
quadrature sensitivity and unresolved derivative-reference study remain in the
[dated evidence](../local/validation-refoundation/2026-09-15/test-system-audit.md).
They are research observations, not fabricated code failures or release approvals.

The optional physical acceptance runner
`test/manual/fem/validate_onelab_export.jl /new/evidence/directory` checks grounded
and insulated bundled terminals, native conductivity/amplitude edits, and
Julia/detached parity after remeshing. It executes real solves and requires
native Gmsh/GetDP executables; it is excluded from automatic discovery.
