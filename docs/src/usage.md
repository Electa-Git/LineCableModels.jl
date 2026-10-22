# Modeling and results

LineCableModels builds complete cable problems from material, geometry, earth,
and frequency data. A scalar declaration constructs one value. A declaration
containing [`Grid`](@ref) or [`Gridspace`](@ref) returns a finite space of
values.

## Computations

The ordinary computation form is `problem → compute` with an optional explicit
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
then returns a target-bearing formulation space. The following declaration
selects both constant earth properties and Portela's frequency-dependent
material law for comparison:

```julia
formulations = Formulation(
    earth_properties = Grid((:constant, :portela1999)),
)

run = compute(
    ParametricProblem(problem_space),
    Combinatorial(formulations),
)
```

The computation contains `length(problem_space) * length(formulations)`
results. Each problem point is materialized once and is evaluated against all
resolved formulations. Use `combine=:product` or `combine=:zip` on the
formulation constructor only to compose fields inside that formulation. The
outer problem-formulation relation is always Cartesian. Formula-owned numerical
controls also vary as complete selections. The following declaration selects
two convergence tolerances for the modal iteration:

```julia
modal_formulations = ModalAnalysisFormulation(Grid((
    formula(:default; options=(iteration=(convergence=1e-4,),)),
    formula(:default; options=(iteration=(convergence=1e-8,),)),
)))
```

Execution settings use the `options` keyword:

```julia
compute(problem, formulation; options=(trace=true,))
```

Line parameters are always per meter, with Z in Ω/m and Y in S/m. The whole line's
matrices are the per-meter ones times the line length stored in the result:

```julia
parameters = compute(problem, formulation)
Z(parameters) .* line_length(parameters)   # Ω
Y(parameters) .* line_length(parameters)   # S
```

Cable constants use a separate earth-free workflow:

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
contain one row per independent concentric assembly. `only(constants)` returns
the scalar row of a conventional single-core coaxial cable. [`LineParameters`](@ref)
stores frequency-dependent Z/Y matrices per meter and records their physical domain
and the line length. Scientific extraction goes through `observe` or `observables`.

[`observe`](@ref) reads one scientific quantity:

```julia
observe(parameters, frequencies)
observe(parameters, R, 1, 1, Colon())
observe(parameters, Z, angle, 1, 1, Colon())
@observe parameters Z[1, 2, :]
```

The [`@observe`](@ref) macro constructs detached requests or immediately
extracts indexed values.

[`observables`](@ref) constructs detached values for tables and plots. Requests
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
label(display_unit(R))               # "Ω/km"
```

## Parameter spaces and uncertainty

`ParametricResult`, `LinearErrorResult`, and `MonteCarloResult` are finite
one-dimensional collections. Indexing and iteration return stored core
results. `first`, `only`, `collect`, `map`, and `zip` retain their ordinary
Julia meanings. A `ParametricResult` also holds the gridpoints of its values and
permits two-axis lookup. Each value records the formulation that computed it:

```julia
run.gridpoints
selected = run[problem_index, formulation_index]
formula_id(details(selected).data.formulation.formulas.earth_impedance)
```

Linear storage is column-major in `(problem, formulation)` coordinates, with the problem index varying fastest. Scalar formulations produce a singleton
formulation axis and preserve the established scalar numerical computation.

When one point of a study fails, the space records it and continues, and the study keeps
the points that it computed. The build or the resolve of a problem, or the computation of a
point, fails the point when it throws an exception that `on_error(study, exception)` returns
for the study. Cables that overlap (`OverlapError`) and an input outside the domain of a formula
(`DomainError`) are infeasible designs, and a study records them. Any other exception stops
the study with its backtrace. A `ParametricResult` and a `LinearErrorResult` hold the values
of the computed points, and `run.failures` holds the gridpoint identity and the exception of
each failed point, in the same order. Indexing, iteration and `length` cover the computed
points, and two-axis lookup of a failed point throws an error that shows its exception. A
study with no computed point throws. A user records an exception type of their own with one
method, `LineCableModels.on_error(::Combinatorial, exception::MyError) = exception`.

Use the product accessors for uncertainty computations:

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
tuple, such as `MonteCarlo(formulation; options=(trials=1000, seed=42))`.
Existing keyword shorthand such as `MonteCarlo(formulation; trials=1000, seed=42)`
uses the same validation. Normalized settings are stored in
`formulation.options`. An options type or cast is unnecessary.

A draw whose build throws `OverlapError`, or whose build or computation throws
`DomainError`, is infeasible. Monte Carlo discards the infeasible draw and draws again,
so that draw enters neither samples nor statistics, and nothing records it. After
`max_failures` infeasible draws in a row for one trial, the point fails with a `DomainError`
whose message gives that count and the last exception. Monte Carlo builds the nominal
problem of every point, from the means of its uncertain arguments, before the first draw,
and a nominal problem that overlaps fails its point. The study records a failed point in
`sampled.failures`, as a `ParametricResult` does, and continues. Any other exception stops
the study. The resulting distribution is conditional on feasible draws.

When `trials=nothing`, Monte Carlo uses a simultaneous
Dvoretzky-Kiefer-Wolfowitz bound. For `M` scalar marginals and confidence
`1-α`, the trial count is
`ceil(log(2M/α) / (2*cdf_tol^2))`. `cdf_tol` bounds empirical-CDF deviation.
It does not bound the mean error or the joint distribution.

## Tables and reports

Pass a completed result directly to [`report`](@ref). The positional selection
takes the same scientific requests as plotting. Omit it to use the result's
available quantities and default display units:

```julia
constants_report = report(constants)
phase_report = report(parameters, (R, L, G, C); length_unit=:kilo)
first_term_report = report(parameters, @observe(R[1, 1, 1:12]))
resistance = phase_report[R]
```

Each quantity has its own DataFrame: full line matrices have one frequency row
per sample and columns for every ordered coefficient, including both
off-diagonals. Cable constants have one operating-frequency row and named
assembly columns. Displaying the report shows these labelled tables in text or
HTML. Collections retain separate tables
for each gridpoint in its original order.

`phase_report[R]` returns the resistance DataFrame already produced by the
report, including any requested coefficient or sample subset. For a collection,
`study_report[R]` returns an ordered vector of DataFrames and `study_report[2, R]`
returns the second result's table. The index is a position in the ordered
observed-result collection, not its original gridpoint identifier. A one-result
collection still returns a vector.
Use `copy(phase_report[R])` when edits should leave the report unchanged.

An indexed request must match the reported selection. Absent or ambiguous
products raise an error.
Select different scientific contents with `report(source, selection)`.

For statistical products, select the statistic explicitly:

```julia
using Statistics
statistics_report = report(sampled, ((statistics, R, mean), (statistics, R, std)))
mean_resistance = statistics_report[1, (statistics, R, mean)]
```

Explicit `ObservedResult` snapshots are useful when retaining selected scientific
products independently of their source. Reporting a snapshot preserves its
recorded units. Supplied compatible unit options re-express the retained values.
Clipping decisions and sample retention are fixed when the observation is
created. Use the report's quantity tables to access these values as DataFrames.

Loading XLSX activates one workbook per gridpoint and quantity. Each destination is checked before writing. Replacing files requires `overwrite=true`:

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
sample means and standard deviations. `uncertain` returns those stored cores.
It does not recover joint output correlations. Loading
Distributions enables supported univariate distributions as Monte Carlo
samplers and `pdf`/`cdf` evaluation for [`HistogramDensity`](@ref).

Repeated use of the same uncertain argument during direct propagation retains
its covariance. Distinct uncertain arguments remain independent.
