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
| `julia --project=test/core test/runtests.jl tag:core_only` | Optional-extension boundaries with core dependencies only. |
| `LINECABLEMODELS_TEST_PLOTTING=true julia --project=test/visual --compiled-modules=no test/runtests.jl tag:visual` | Existing rendering, plot lifecycle, layout and output behaviour checks; no new calibration. |
| `julia --project=test test/runtests.jl tag:pscad_native` | Real PSCAD acceptance. Set `LINECABLEMODELS_PSCAD_CONFIG` to your private TOML filename; requires the configured live station. |
| `julia --project=test test/runtests.jl tag:aqua` | Aqua in a fresh Julia process. |
| `julia --project=docs docs/doctest.jl` and `julia --project=docs docs/make.jl` | Docstrings and documentation build; instantiate with `docs/instantiate.jl`. |
| `julia --project=ext/LineCableModelsPSCADExt/remote test/integration/pscad/remote_runner_boundary.jl` | Remote-runner protocol against a mocked Python automation boundary; does not execute PSCAD or validate solver numerics. |

Ordinary exclusions are `quality`, `aqua`, `visual`, `core_only` and `pscad_native`.

The FEM extension tests cover native parameters, geometry, mesh inspection,
export, solver orchestration and result import. `fem_native_mesh.jl` evaluates
both native parsers and generates meshes without field solves. Independent
manufactured Maxwell fixtures check the field equations. Numerical integration
tests remain separate from the manual FEM/reference comparisons in
`test/manual/calculations/fem_validation/`.

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
Then merge and enforce the gate:

```sh
julia --project=test/coverage test/coverage.jl check
```

This inventories every production Julia file under `src/` and `ext/`, writes
`lcov.info`, and fails below 95%. The check removes traces afterward, including on failure;
retain the report. Missing executions and failing tests remain visible even if
the line percentage passes.

## Verification limits

Software tests check implemented behaviour checks. FEM/reference comparisons and mesh
convergence need the explicit manual tool and numerical evidence. A passed
software suite or coverage gate does not qualify a frequency range or material.
The PSCAD acceptance items need the configured private station. Graphics checks
need their backend environment and a display where the backend needs one.
