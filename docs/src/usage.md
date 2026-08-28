# Modelling and results

LineCableModels builds complete cable problems from material, geometry, earth,
and frequency data. A scalar declaration constructs one value. A declaration
containing [`Grid`](@ref) or [`Gridspace`](@ref) returns a finite space of
values.

## Calculations

The ordinary calculation form is `problem → compute` with an optional explicit
formulation:

```julia
parameters = compute(problem)
parameters = compute(problem, Formulation())

all_points = compute(
    ParametricProblem(problem_space),
    Combinatorial(Formulation()),
)

linear = compute(
    ParametricProblem(problem_space),
    LinearError(Formulation()),
)

sampled = compute(
    ParametricProblem(problem_space),
    MonteCarlo(Formulation(); trials=1000, seed=42),
)
```

[`Formulation`](@ref) selects the physical and numerical methods used for one
problem. `Combinatorial` evaluates every selected point. `LinearError` applies
direct linear uncertainty propagation. `MonteCarlo` samples independent
realisations within each selected point.

Execution settings use the `options` keyword:

```julia
compute(problem, formulation; options=(output_basis=:total,))
```

`output_basis=:total` scales both impedance and admittance by the system
length. `CableConstants(design)` uses a 1 m ordinary line problem and returns
per-length values.

## Completed results

[`CableConstants`](@ref) stores R/L/C values per metre. [`LineParameters`](@ref)
stores frequency-dependent Z/Y matrices and records their physical domain and
`:pul` or `:total` basis. Scientific extraction goes through `observe` or
`observables`.

[`observe`](@ref) reads one scientific quantity:

```julia
observe(parameters, frequencies)
observe(parameters, R, 1, 1, Colon())
observe(parameters, Z, angle, 1, 1, Colon())
@observe parameters Z[1, 2, :]
```

The [`@observe`](@ref) macro constructs detached requests or immediately
extracts indexed values.

[`observables`](@ref) prepares detached values for tables and plots. Requests
are positional:

```julia
published = observables(
    parameters,
    (
        @observe(R[1, 1, :]),
        @observe((Z, angle)[1, 1, :]),
    ),
)
```

Each returned item has `values`, `quantity`, and `unit`. Labels, symbols, and
display units come from the quantity metadata:

```julia
label(R)                             # "Series resistance"
symbol(Z, angle)                     # "∠Z"
label(display_unit(R, :pul))         # "Ω/km"
```

## Parameter spaces and uncertainty

`ParametricResult`, `LinearErrorResult`, and `MonteCarloResult` are finite
one-dimensional collections in Gridspace traversal order. Indexing and
iteration return stored core results. `first`, `only`, `collect`, `map`, and
`zip` retain their ordinary Julia meanings.

Use the product accessors for uncertainty calculations:

```julia
result(sampled)
statistics(sampled)
samples(sampled)
histograms(sampled)
uncertain(sampled)
```

Monte Carlo settings and resolved point data are available through
`root_seed`, `point_seed`, `trial_count`, `confidence`, `cdf_tolerance`, and
`sampling_distribution`.

The default error mode propagates every exception. Conditional rejection of
unsupported realisations is explicit:

```julia
MonteCarlo(
    formulation;
    trials=1000,
    options=(
        retain_details=true,
        on_error=:retry,
        max_failures=100,
    ),
)
```

Only `DomainError` is retryable. Rejected values do not enter samples or
statistics. `details(result).failure_summary` reports attempts, accepted and
failed counts, acceptance rate, and failure counts by error type and stage.
The resulting distribution is conditional on successful problem construction
and calculation.

When `trials=nothing`, Monte Carlo uses a simultaneous
Dvoretzky–Kiefer–Wolfowitz bound. For `M` scalar marginals and confidence
`1-α`, the trial count is
`ceil(log(2M/α) / (2*cdf_tol^2))`. `cdf_tol` bounds empirical-CDF deviation;
it does not bound the mean error or the joint distribution.

## Tables and reports

DataFrames constructors return one wide table. Coordinates identify rows and
each requested quantity occupies one column:

```julia
using DataFrames

table = DataFrame(parameters, (@observe(R[:, :, :]), @observe(L[:, :, :])))
summary = DataFrame(sampled)
```

[`report`](@ref) creates an in-memory table or a written report from explicit
requests. Loading XLSX activates workbook output:

```julia
using XLSX

artifact = report(
    XLSXReportDefinition(file_name="line_parameters.xlsx"),
    parameters,
)
artifact.output
```

## Plots

Load one Makie package before calling `plot` or `preview`:

```julia
using CairoMakie

plots = plot(parameters)
preview(cable)
```

Monte Carlo results use native Makie verbs and an explicit marginal:

```julia
Makie.hist(sampled, R; bins=20)
Makie.ecdfplot(sampled, R)
Makie.qqplot(sampled, R; qqline=:identity)

Makie.hist(sampled, @observe(L[1, 1, 3]); bins=20)
```

Each non-mutating managed plotting call returns [`UIPlot`](@ref). Use
[`export_svg`](@ref) to save the current view.

## Optional uncertainty packages

The core `Grid`/`Gridspace` grammar does not load Measurements.jl or
Distributions.jl. Loading Measurements enables `LinearError`. Loading
Distributions enables supported univariate distributions as Monte Carlo
samplers and `pdf`/`cdf` evaluation for [`HistogramDensity`](@ref).

Repeated use of the same uncertain argument during direct propagation retains
its covariance. Distinct uncertain arguments remain independent.
