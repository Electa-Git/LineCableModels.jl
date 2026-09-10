# Test suite

Tests assert public results and invariants rather than copying source functions
or tutorial scripts. Files use lowercase snake case and mirror the source
module they exercise.

Run the default unit, integration, and non-graphical extension suite with:

```julia
using Pkg
Pkg.test()
```

The package and `test/Project.toml` form one Julia 1.12 workspace. Add ordinary
test-only dependencies to that test project. The visual, core-only, gauntlet,
and coverage projects remain isolated because they verify different dependency
boundaries.

Pass selectors through `test_args`. A plain selector matches a file path or test-item
name. A `tag:` selector matches a test tag:

```julia
Pkg.test(test_args = ["tag:unit"])
Pkg.test(test_args = ["tag:integration"])
Pkg.test(test_args = ["Engine / solver"])
```

The package CLI has independent process-boundary tests requiring only Python's
standard library:

```sh
python3 -m unittest discover -s test/cli -v
```

They cover symlink installation and invocation, application selection, arguments,
working directories, exit codes and termination. The quality CI job runs them
without Julia packages or application services.

The supported tags are `unit`, `integration`, `extension`, `fem_numerical`, `visual`, `quality`, `gauntlet`, and `gauntlet_toolkit`. Visual, quality, `core_only`, and both gauntlet tags are excluded from the default run and execute in dedicated environments. See
[`gauntlet/README.md`](../gauntlet/README.md) for explicit campaign, comparison and recovery commands.

The `fem_numerical` items exercise deterministic multi-frequency FEM solves and
the frozen Python reference matrices. They are excluded from the ordinary test
run and execute explicitly with the package's hash-pinned GetDP 3.5.0 complex
artifact. To run them locally:

```sh
julia --project=test \
  -e 'push!(ARGS, "tag:fem_numerical"); include("test/runtests.jl")'
```

Set `LINECABLEMODELS_GETDP=/absolute/path/to/getdp` only to exercise an external
solver override.

The full gauntlet remains a manual workflow, separate from these deterministic
FEM regressions. CI must not launch live PSCAD/FEM gauntlet campaigns or promote
their output to references. The read-only [numerical-reference gate](numerical/README.md)
replays stored problem and formulation declarations against explicitly reviewed,
pinned gauntlet arrays, with element-wise RMS tolerances and scalar inference
checks. Its approval manifest and bindings remain empty; it is not enabled in CI.

Instantiate the gauntlet environment and run every tagged case through the dedicated TestItemRunner entry point:

```sh
julia --project=gauntlet -e 'using Pkg; Pkg.instantiate()'
julia --project=gauntlet \
  test/gauntlet/runtests.jl
```

`test/gauntlet/runtests.jl` selects both `gauntlet` and `gauntlet_toolkit` tests.
These check declarations, execution, persistence and comparison without native
PSCAD. Reusable declarations live in `gauntlet/benchmarks/`; indexed physical
models live in `gauntlet/cases/`. Live campaigns are explicit CLI actions.

During development, select the owned UQ or external PSCAD family directly:

```sh
julia --project=gauntlet --startup-file=no -e \
  'using TestItemRunner; TestItemRunner.run_tests(joinpath(pwd(), "test"); filter=ti -> :uq in ti.tags, verbose=true)'
julia --project=gauntlet --startup-file=no -e \
  'using TestItemRunner; TestItemRunner.run_tests(joinpath(pwd(), "test"); filter=ti -> :pscad in ti.tags, verbose=true)'
```

Run the reusable gauntlet toolkit checks separately with:

```sh
julia --project=gauntlet --startup-file=no \
  -e 'push!(ARGS, "tag:gauntlet_toolkit"); include("test/runtests.jl")'
```

The PSCAD worker uses its existing PythonCall environment for local protocol tests.
The automation double checks settings, file collection and cleanup; it does not
simulate electromagnetic results or require a PSCAD installation:

```sh
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE=python3 \
  julia --project=ext/LineCableModelsPSCADExt/remote -e 'using Pkg; Pkg.instantiate()'
JULIA_CONDAPKG_BACKEND=Null JULIA_PYTHONCALL_EXE=python3 \
  julia --project=ext/LineCableModelsPSCADExt/remote test/integration/pscad/worker.jl
```

The deterministic Cairo suite has its own environment and may be run headlessly with:

```sh
julia --project=test/visual -e 'using Pkg; Pkg.instantiate()'
LINECABLEMODELS_TEST_PLOTTING=true julia --project=test/visual \
  -e 'push!(ARGS, "tag:visual"); include("test/runtests.jl")'
```

Exercise optional-dependency fallbacks without loading any weak dependency:

```sh
julia --project=test/core -e 'using Pkg; Pkg.instantiate()'
julia --project=test/core \
  -e 'push!(ARGS, "tag:core_only"); include("test/runtests.jl")'
```

Fixture factories live in `support/fixtures.jl` and must return fresh mutable objects.
Static input data belongs in `fixtures/data`, independently sourced numerical values in
`fixtures/reference`, and current rendering baselines in `fixtures/golden`. Ordinary
tests write only to system temporary directories. Explicit live Gauntlet runs preserve
projects and diagnostics in the declared campaign and backend work directories.

Numerical tests use scale- and precision-aware helpers from `support/numerical.jl`.
There is no suite-wide absolute tolerance. Expected values must come from analytical
identities, independent references, residuals, or other observable invariants—not from
reimplementing the function under test.

Regression tests preserve current scientific and architectural behavior, not every
spelling used during unreleased development:

- Keep exact checks for identity, units, ordering, and repeated scalar/batched
  calculations in the same runtime.
- Compare archived floating-point references component by component at their own
  scale. For ill-conditioned equivalent permeability, preserve the physical GMR
  rather than the last bits of the intermediate coefficient. Do not refresh the
  frozen reference merely because a dependency changes final rounding.
- Treat recovered formula inventories as required subsets. Inspect every currently
  registered formula for route ownership and documentation without prohibiting new
  registrations.
- Keep absence tests for deliberately retired abstractions. Exercise orchestration
  counts/order at runtime; private helper names and import formatting are not APIs.
- Mesh every supported bounded formation and check elements on every material
  surface. A nonempty overall mesh or physical-tag list is insufficient.

The enforced coverage ratio includes all production code under `src/` and
`ext/`. The LCOV report also publishes reusable gauntlet helper coverage when traces
exist, while excluding manually authored files under `gauntlet/cases/`. Clean stale
traces before a coverage run, merge traces from the ordinary, core-only, and visual
environments, and enforce the source-amended 95% gate afterward:

```sh
julia --project=test/coverage -e 'using Pkg; Pkg.instantiate()'
julia --project=test/coverage test/coverage.jl clean
julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'
julia --project=test/core --code-coverage=@. \
  -e 'push!(ARGS, "tag:core_only"); include("test/runtests.jl")'
LINECABLEMODELS_TEST_PLOTTING=true julia --project=test/visual \
  --compiled-modules=no --code-coverage=@. \
  -e 'append!(ARGS, ["tag:visual", "loaded extension activation"]); include("test/runtests.jl")'
julia --project=test/coverage test/coverage.jl check
```

CI additionally collects the local PSCAD worker protocol traces and checks backend selection in isolated GLMakie and
WGLMakie environments (GLMakie runs under Xvfb). The deterministic FEM job uploads
its production traces with a `.fem.cov` suffix to avoid cross-runner process-ID
collisions. The solver-free toolkit job uploads `.toolkit.cov` traces as well;
its reusable helper coverage is published, but only `src/` and `ext/` contribute
to the production threshold. Codecov's project and changed-line checks use these
same production paths and retain their 95% targets. Gauntlet files remain visible
in the published report. Neither job runs the manual Gauntlet campaign.
The coverage job merges both artifacts before the same single check. The cleaner
removes Julia-native, imported FEM, imported toolkit and documentation traces after
all instrumented workers exit, even when the coverage check fails. The checker amends
coverage from source, rejects any missing `src/` or `ext/` Julia file, writes
`lcov.info`, and fails below 95% aggregate line coverage.

Documentation, Aqua, and golden regeneration are separate checks and must not
be used to satisfy the production coverage threshold.
