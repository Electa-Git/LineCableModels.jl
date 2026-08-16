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

The supported tags are `unit`, `integration`, `extension`, `visual`, and `quality`.
Visual, quality, and `core_only` boundary items are intentionally excluded from the
default run and execute in dedicated CI environments.

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
`fixtures/reference`, and current rendering baselines in `fixtures/golden`. Tests must
write only to system temporary directories.

Numerical tests use scale- and precision-aware helpers from `support/numerical.jl`.
There is no suite-wide absolute tolerance. Expected values must come from analytical
identities, independent references, residuals, or other observable invariants—not from
reimplementing the function under test.

Coverage is collected only for tracked production code under `src/` and `ext/`. Clean
stale traces before a coverage run, merge traces from the ordinary, core-only, and
visual environments, and enforce the source-amended 95% gate afterward:

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
