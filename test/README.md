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

The supported tags are `unit`, `integration`, `extension`, `fem_numerical`, `visual`, `quality`, `gauntlet`, and `gauntlet_toolkit`. Visual, quality, `core_only`, and both gauntlet tags are excluded from the default run and execute in dedicated environments. See
[`gauntlet/README.md`](gauntlet/README.md) for the explicit snapshot, live, and record
commands.

The `fem_numerical` items exercise deterministic multi-frequency FEM solves and
the frozen Python reference matrices. Without GetDP they report explicit skips,
not passing assertions. CI runs them in a separate job with a checksum-pinned
GetDP 3.5.0 complex solver and fails setup if that executable is unavailable.
To run those checks locally with an installed solver:

```sh
LINECABLEMODELS_GETDP=/absolute/path/to/getdp julia --project=test \
  -e 'push!(ARGS, "tag:fem_numerical"); include("test/runtests.jl")'
```

The full gauntlet remains a manual workflow, separate from these deterministic
FEM regressions. CI must not launch live PSCAD/FEM gauntlet campaigns or promote
their output to references. A future artifact-backed gate must consume explicitly
validated, pinned gauntlet artifacts; that binding is not configured yet.

Instantiate the gauntlet environment and run every tagged case through the dedicated TestItemRunner entry point:

```sh
julia --project=test/gauntlet -e 'using Pkg; Pkg.instantiate()'
LINECABLEMODELS_GAUNTLET_MODE=snapshot julia --project=test/gauntlet \
  test/gauntlet/runtests.jl
```

`test/gauntlet/runtests.jl` uses `@run_package_tests(filter=ti -> :gauntlet in ti.tags, verbose=true)`. TestItemRunner discovers every benchmark file below `gauntlet/benchmarks/`; the indexed, backend-neutral physical models live separately below `gauntlet/cases/`. The command does not maintain a second benchmark list. Each benchmark file contains exactly one benchmark, an invariant enforced by the toolkit suite.

During development, select the owned UQ or external PSCAD family directly:

```sh
julia --project=test/gauntlet --startup-file=no -e \
  'using TestItemRunner; TestItemRunner.run_tests(joinpath(pwd(), "test"); filter=ti -> :uq in ti.tags, verbose=true)'
julia --project=test/gauntlet --startup-file=no -e \
  'using TestItemRunner; TestItemRunner.run_tests(joinpath(pwd(), "test"); filter=ti -> :pscad in ti.tags, verbose=true)'
```

Run the reusable gauntlet toolkit checks separately with:

```sh
julia --project=test/gauntlet --startup-file=no \
  -e 'push!(ARGS, "tag:gauntlet_toolkit"); include("test/runtests.jl")'
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
tests write only to system temporary directories. Explicit live gauntlet runs preserve
their disposable projects and diagnostics below `test/gauntlet/benchmarks/.work/`.

Numerical tests use scale- and precision-aware helpers from `support/numerical.jl`.
There is no suite-wide absolute tolerance. Expected values must come from analytical
identities, independent references, residuals, or other observable invariants—not from
reimplementing the function under test.

The enforced coverage ratio includes only tracked production code under `src/` and
`ext/`. The LCOV report also publishes reusable gauntlet helper coverage when traces
exist, while excluding manually authored files under `test/gauntlet/cases/`. Clean stale
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

CI additionally checks backend selection in isolated GLMakie and
WGLMakie environments before the same single merge/check step. The checker amends
coverage from source, rejects any missing `src/` or `ext/` Julia file, writes
`lcov.info`, and fails below 95% aggregate line coverage.

Documentation, Aqua, and golden regeneration are separate checks and must not
be used to satisfy the production coverage threshold.
