# Modeling and results

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
realizations within each selected point.

Any final formulation slot can be an explicit finite source. The constructor
then returns a target-bearing formulation space. For example, compare constant
earth properties with Portela's frequency-dependent material law:

```julia
formulations = Formulation(
    earth_properties = Grid((:constant, :portela1999)),
)

run = compute(
    ParametricProblem(problem_space),
    Combinatorial(formulations),
)
```

The calculation contains `length(problem_space) * length(formulations)`
results. Each problem point is materialized once and is evaluated against all
resolved formulations. Use `combine=:product` or `combine=:zip` on the
formulation constructor only to compose fields inside that formulation; the
outer problem/formulation relation is always Cartesian. Formula-owned numerical
controls also vary as complete selections; for example, modal iteration convergence:

```julia
modal_formulations = ModalTransformationFormulation(Grid((
    formula(:default; options=(iteration=(convergence=1e-4,),)),
    formula(:default; options=(iteration=(convergence=1e-8,),)),
)))
```

Execution settings use the `options` keyword:

```julia
compute(problem, formulation; options=(output_basis=:total,))
```

`output_basis=:total` scales both impedance and admittance by the system
length. Cable constants use a separate earth-free workflow:

```julia
constants_problem = CableConstantsProblem(
    design;
    temperature = 20.0,
    frequency = 50.0,
)
constants = compute(constants_problem, CableConstantsFormulation())

# Admitted convenience for the same operation
constants = CableConstants(design; temperature = 20.0, frequency = 50.0)
```

The cable-constant workflow admits 50 Hz or 60 Hz and has no earth model,
placement, propagation constant, transposition, or bundle option.

## Completed results

[`CableConstants`](@ref) stores R/L/C/G values per meter. Its aligned vectors
contain one row per independent concentric assembly; `only(constants)` returns
the scalar row of a conventional single-core coaxial cable. [`LineParameters`](@ref)
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
one-dimensional collections. Indexing and iteration return stored core
results. `first`, `only`, `collect`, `map`, and `zip` retain their ordinary
Julia meanings. A `ParametricResult` also retains its resolved axes and permits
two-axis lookup:

```julia
run.axes.problems
run.axes.formulations
selected = run[problem_index, formulation_index]
formula_id(run.axes.formulations[formulation_index].methods.earth_impedance)
```

Linear storage is column-major in `(problem, formulation)` coordinates, so the
problem index varies fastest. Scalar formulations produce a singleton
formulation axis and preserve the established scalar numerical calculation.

Use the product accessors for uncertainty calculations:

```julia
collect(sampled)
statistics(sampled)
samples(sampled)
histograms(sampled)
uncertain(sampled)
```

Monte Carlo settings and resolved point data are available through
`root_seed`, `point_seed`, `trial_count`, `confidence`, `cdf_tolerance`, and
`sampling_distribution`.

Every Monte Carlo execution control may be supplied in an ordinary named
tuple, for example `MonteCarlo(formulation; options=(trials=1000, seed=42))`.
Existing keyword shorthand such as `MonteCarlo(formulation; trials=1000, seed=42)`
uses the same validation. Normalized settings are stored in
`formulation.options`; no options type or cast is required.

The default error mode propagates every exception. Conditional rejection of
unsupported realizations is explicit:

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
statistics. `details(result).data.failure_summary` reports attempts, accepted and
failed counts, acceptance rate, and failure counts by error type and stage.
The resulting distribution is conditional on successful problem construction
and calculation.

When `trials=nothing`, Monte Carlo uses a simultaneous
Dvoretzky–Kiefer–Wolfowitz bound. For `M` scalar marginals and confidence
`1-α`, the trial count is
`ceil(log(2M/α) / (2*cdf_tol^2))`. `cdf_tol` bounds empirical-CDF deviation;
it does not bound the mean error or the joint distribution.

## Tables and reports

Pass a completed result directly to [`report`](@ref). `values` accepts the same
scientific requests that plotting names `ydata`. Omit it to use the result's
available quantities and default display units:

```julia
constants_report = report(constants)
phase_report = report(parameters; values=(R, L, G, C), length_unit=:kilo)
first_term_report = report(parameters; values=@observe(R[1, 1, 1:12]))
resistance = phase_report[R]
```

Each quantity has its own DataFrame: full line matrices have one frequency row
per sample and columns for every ordered coefficient, including both
off-diagonals. Cable constants have one operating-frequency row and named
assembly columns. Displaying the report shows these labelled tables in text or
HTML. Collections retain separate tables
for each gridpoint in its original order. Default reporting creates no figure
and writes no files.

`phase_report[R]` returns the resistance DataFrame already produced by the
report, including any requested coefficient or sample subset. For a collection,
`study_report[R]` returns an ordered vector of DataFrames and `study_report[2, R]`
returns the second result's table. The index is a position in the ordered
observed-result collection, not its original gridpoint identifier. A one-result
collection still returns a vector.
Use `copy(phase_report[R])` when edits should leave the report unchanged.

Lookup does not perform further selection or calculation. An indexed request
must match the reported selection; absent or ambiguous products raise an error.
Select different scientific contents with `report(...; values=...)`.

For statistical products, select the statistic explicitly:

```julia
using Statistics
statistics_report = report(sampled; values=((statistics, R, mean), (statistics, R, std)))
mean_resistance = statistics_report[1, (statistics, R, mean)]
```

Explicit `ObservedResult` snapshots are useful when retaining selected scientific
products independently of their source. Reporting a snapshot preserves its
recorded units; supplied compatible unit options re-express the retained values.
It cannot revise prior clipping or recover discarded samples. `DataFrame(observed)`
is not an aggregate conversion; inspect the report's quantity tables instead.

Loading XLSX activates one workbook per gridpoint and quantity. All
destinations are checked before writing; replacing files requires `overwrite=true`:

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
Distributions.jl. Load Measurements before `LinearError` or `MonteCarlo`
computation. MC constructs and stores marginal Measurements from accepted raw
sample means and standard deviations; `uncertain` returns those stored cores.
It does not recover joint output correlations. Loading
Distributions enables supported univariate distributions as Monte Carlo
samplers and `pdf`/`cdf` evaluation for [`HistogramDensity`](@ref).

Repeated use of the same uncertain argument during direct propagation retains
its covariance. Distinct uncertain arguments remain independent.
