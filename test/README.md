# Test suite

The suite is organized around observable contracts rather than source functions or
tutorial scripts. Files use lowercase snake case and mirror the production subsystem
they exercise.

Run the default unit, integration, and non-graphical extension suite with:

```julia
using Pkg
Pkg.test()
```

Pass selectors through `test_args`. A plain selector matches a file path or test-item
name; a `tag:` selector matches a test tag:

```julia
Pkg.test(test_args = ["tag:unit"])
Pkg.test(test_args = ["tag:integration"])
Pkg.test(test_args = ["Engine / solver"])
```

The supported tags are `unit`, `integration`, `extension`, `visual`, `quality`, and
`gauntlet`. Visual, quality, `core_only`, and gauntlet items are intentionally excluded
from the default run and execute in dedicated environments. See
[`gauntlet/README.md`](gauntlet/README.md) for the explicit snapshot, live, and record
commands.

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
their disposable projects and diagnostics below `test/gauntlet/cases/.work/`.

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
LINECABLEMODELS_TEST_PLOTTING=true julia --project=test/visual --code-coverage=@. \
  -e 'append!(ARGS, ["tag:visual", "loaded extension activation"]); include("test/runtests.jl")'
julia --project=test/coverage test/coverage.jl check
```

CI additionally collects the backend-selection contract in isolated GLMakie and
WGLMakie environments before the same single merge/check step. The checker amends
coverage from source, rejects any missing `src/` or `ext/` Julia file, writes
`lcov.info`, and fails below 95% aggregate line coverage.

Documentation, Aqua, and golden regeneration are separate correctness gates and must not
be used to satisfy the production coverage threshold.
