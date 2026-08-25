# API reference

## Modelling and execution grammar

The top-level namespace exposes the declarative modelling API. `Grid` is the
only parameter-variation syntax, and `Gridspace{Target}` is a lazy space of
complete targets. Scalar `Material` keywords return one materialised value,
while varying inputs return a `Gridspace`. `CableBuilder`, `Earth`, and
`SystemBuilder` likewise return eager domain values for scalar input and a
`Gridspace` only when an input is explicitly varied.

Ordinary collections remain ordinary constructor values. Wrap a collection in
`Grid` only when it is intended to vary. Reusing the same `Grid` object has no
hidden effect: `:product` still forms the Cartesian product and `:zip` still
pairs rows.

```julia
earth = Earth(
    rho=Grid((10.0, 100.0, 1000.0)),
    eps_r=Grid((100.0, 10.0, 5.0));
    combine=:zip,
)
```

Composition is local to each Gridspace node. The example above is already a
space with exactly three earth points. The space can still participate in a
Cartesian product at a parent problem node.

See the [Gridspace manual](gridspace.md) for composition, eager construction,
uncertainty realisation, and extension methods.

All numerical work follows `problem → formulation → compute`:

Here, `problem_space` is the `Gridspace` returned by `SystemBuilder`.

```julia
compute(first(problem_space), Formulation())

compute(ParametricProblem(problem_space), Combinatorial(Formulation()))
compute(ParametricProblem(problem_space), LinearError(Formulation()))
compute(
    ParametricProblem(problem_space),
    MonteCarlo(Formulation(); trials=1000, cdf_tol=0.02),
)
```

The higher-order formulation is always explicit. `Combinatorial` evaluates
every selected point, `LinearError` selects direct linear uncertainty propagation,
and `MonteCarlo` samples realisations within each selected outer point. Every
cardinality, including one, returns the corresponding composite result family.

`Formulation` owns physics and numerical-method choices. A computation owner
validates its execution tuple through [`computation_options`](@ref). A
composite calculation routes separate tuples when it invokes different
backends. See [Computational engine](@ref) for the option API and extension
methods. A
materialised line system can request total rather than per-length Z/Y matrices
without altering its formulation:

```julia
compute(problem, formulation; options=(output_basis=:total,))
```

`output_basis=:total` scales both impedance and admittance by the materialised
system length. `CableConstantsProblem` has no line length and therefore accepts
only the default `:pul` basis.

## Results

`CableConstants` stores R/L/C values per metre. `LineParameters` stores
frequency-dependent Z/Y matrices with their domain and `:pul` or
`:total` basis.

[`observe`](@ref) reads one native numerical meaning from a completed result.
Selectors are the existing function objects:

```julia
observe(parameters, frequencies)
observe(parameters, R, 1, 1, Colon())
observe(parameters, Z, angle, 1, 1, Colon())
@observe parameters Z[1, 2, :]
```

The laconic `Z`, `Y`, `R`, `X`, `L`, `G`, `B`, and `C` accessors delegate to
the same methods. `observe` performs no display conversion and returns no
metadata wrapper. Direct `parameters.Z`, `parameters.Y`, indexing, and views
remain supported for ordinary numerical work. `@observe` is syntax for one
three-index `observe` call; it does not consult the publication declaration.

[`observables`](@ref) is the detached presentation boundary. It requires an
explicit named tuple of requests:

```julia
published = observables(
    parameters,
    (
        frequency = (frequencies, Colon()),
        resistance = (R, 1, 1, Colon()),
        phase = (Z, angle, 1, 1, Colon()),
    ),
)
```

Each field has exactly `values`, `quantity`, and `unit`. `values` is detached
from the result and expressed in `unit`; `quantity` is its precise scientific
identity. Labels and symbols are derived from `LineCableModels.Units`, not
copied into the payload. `observables(typeof(result))` declares the supported
selector vocabulary. There is no zero-argument result publication method.

The same selectors provide quantity metadata directly:

```julia
using LineCableModels

label(R)                             # "Series resistance"
symbol(R)                            # "R"
label(display_unit(R, :pul))         # "Ω/km"

label(Z, abs)                        # "Series impedance magnitude"
symbol(Z, angle)                     # "∠Z"
label(display_unit(Z, angle, :pul))  # "°"
```

`quantity(R)` returns the fieldless typed identity
`LineCableModels.Units.Quantity{:series_resistance}()`. This type is used by
extensions and publication payloads; ordinary result access continues to use
`R`, `Z`, and the other scientific selector functions.

A complete deterministic traversal returns `ParametricResult{T}`. Direct
propagation returns `LinearErrorResult{T}`, and conditional sampling returns
`MonteCarloResult{T}`. In every family, `T` is the primitive `CableConstants`
or `LineParameters` result rather than another composite result.

Use `result`, `statistics`, `samples`, `histograms`, and `uncertain_value` to
inspect scientific products. For a Monte Carlo result, `root_seed`,
`point_seed`, `trial_count`, `confidence`, `cdf_tolerance`, and
`sampling_distribution` read the calculation settings and the values resolved
for each Gridspace point. All three higher-order result families also store a
concrete [`ComputationDetails`](@ref) named tuple. [`details`](@ref) returns
that tuple. It is `(; )` unless the higher-order formulation was constructed
with `options=(retain_details=true,)` and the primitive computation owner
implements [`computation_details`](@ref).

Parametric and linear results retain one detail record per primitive result
under `details(result).points`. Monte Carlo retains one vector per Gridspace
point and one record per trial under `details(result).trials`. Typed details do
not replace statistics, samples, histograms, the root seed, point seeds, or
trial counts. No completed result stores the temporary Gridspace point or a
copy of traversal internals.

`primitives` and `preprocess` are reserved action generics for explicitly
selected future calculation orderings. LineCableModels defines no methods for
either generic: there is no zero-argument conversion, broad
fallback, or implicit transformation. Unsupported orderings fail through
ordinary Julia dispatch.

[`report`](@ref) builds human-facing tables from explicit publication requests.
Its fixed sequence is `entitle → select → tabulate → illustrate →
encode → write → finish`. Existing `DataFrame` methods for completed
results delegate to this boundary. `DataFrame(result::MonteCarloResult)` renders
stored marginal summaries without repeating the calculation. Cable-constant
results produce one R/L/C table. Line-parameter results produce one R/L/C/G
table for every matrix entry and frequency. The displayed `confidence` and
`cdf_tol` values describe the DKW bound below and are not confidence intervals
for the sample mean.

`report(TableReport(...), source)` returns an in-memory
[`ReportArtifact`](@ref) with `output === nothing`. Human-facing line-parameter
workbooks use [`XLSXReport`](@ref):

```julia
using XLSX

artifact = report(
    XLSXReport(file_name="line_parameters.xlsx"),
    parameters,
)
artifact.output
```

Loading XLSX activates the writer extension. The retained
`export_data(:xlsx, parameters; ...)` convenience call delegates to this report
and returns `artifact.output`. ImportExport contains no separate XLSX workbook
path.

After loading a Makie package, retained samples and histograms can be displayed
through the maintained Monte Carlo recipe:

```julia
using CairoMakie

plot(result, R; mode=:hist, data=:both)
plot(result, R; mode=:pdf)
plot(result, R; mode=:ecdf, data=:both)
plot(result, R; mode=:qq)

# Select one line-parameter matrix entry and frequency:
plot(line_result, L; ijk=(1, 1, 3), mode=:hist)
```

The `:pdf`, `:ecdf`, and `:qq` views use the retained piecewise-constant
`HistogramDensity`. When only samples were retained, the recipe derives the
histogram needed for presentation.

When `MonteCarlo(trials=nothing)` is used, the trial count follows a simultaneous
Dvoretzky–Kiefer–Wolfowitz bound. For `M` scalar marginals and confidence
`1-α`, the implementation selects `ceil(log(2M/α) / (2*cdf_tol^2))` trials.
A union bound therefore controls the largest empirical-CDF deviation among
those marginals. `M` is three for cable constants and four real R/L/G/C
upper-triangular coordinates per frequency for line parameters. `cdf_tol` is
not a solver tolerance and does not bound mean error or the joint distribution.

## Optional uncertainty packages

The core Grid/Gridspace grammar does not load Measurements.jl or
Distributions.jl. Loading Measurements adds direct `LinearError(inner)`
propagation. When one selected uncertain argument is passed once to a callable
builder and used more than once there, all uses preserve that exact variable's
covariance. Distinct uncertain arguments remain independent.
`MonteCarloResult` has no implicit conversion to Measurements values.
Empirical reconstruction is an explicit consumer operation and requires
retained joint samples. Loading Distributions adds compatible univariate
distributions as Monte Carlo samplers and `pdf`/`cdf` evaluation of
`HistogramDensity`. A supplied distribution is standardised to the Grid
descriptor's nominal value and standard deviation, so it must have a finite
mean and a positive finite standard deviation.

## Strict materialised API

The eager cable and system types are not exported at top level. Import them
explicitly when strict construction is required. `Material`,
`MaterialsLibrary`, and `CablesLibrary` are public
because they are shared by both modelling modes:

```julia
using LineCableModels: Material
import LineCableModels.DataModel
import LineCableModels.EarthProps

copper = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
core = DataModel.Tubular(0.0, 10e-3, copper)
```

Materialised radial constructors accept resolved numeric inner and outer radii.
Radius-versus-thickness intent and repetition belong to builder definitions.
Group `add!` accepts a fully materialised part of the same scalar type as its
destination.

All eager cable quantities are evaluated at the common material reference
temperature. The materialised `CableDesign` therefore rejects layers with
different `T0` values. Operating temperature belongs only to
`LineParametersProblem`. `SystemBuilder` exposes that setting and [`compute`](@ref)
applies the resistivity correction locally.

## Static earth and solver inspection

`EarthModel` stores only the physical layer description. Frequencies belong to
`LineParametersProblem`, while frequency-dependent earth-property selection
belongs to the analytical formulation:

```julia
using LineCableModels
import LineCableModels.Engine

formulation = Formulation(
    :analytical;
    earth_properties=Engine.EarthProperties.CPEarth(),
)
```

Ordinary `compute` returns `LineParameters`. Select trace output on the
formulation when completed internal matrices are required:

```julia
trace_formulation = Formulation(
    :analytical;
    earth_properties=Engine.EarthProperties.CPEarth(),
    options=(output=:trace,),
)
trace = compute(problem, trace_formulation)
parameters = trace.result
```

The trace owns `Zin`, `Pin`, `Zg`, `Pg`, `Z`, and `P` tensors. Parameter-only
calculation does not allocate these three-dimensional snapshots.

Component-specific verbosity is configured with the computation-options named
tuple. A missing key uses `default`, and logging remains scoped to the
calculation:

```julia
options = (verbosity=(default=0, LineCableModels=1, NLsolve=0, QuadGK=0),)
compute(problem, formulation; options)
```

## Wire estimates

[`make_stranded`](@ref) and [`make_screened`](@ref) return a [`WireEstimate`](@ref).
When no candidate meets every requested limit, the search returns ranked
candidates and records the unmet constraints:

```julia
estimate = make_stranded(1000.0)
closest = estimate[:match]
fewest_layers = estimate[:layers]
```

The `:match`, `:layers`, `:wires`, and `:diameter` selectors use dispatch and never
rewrite the retained candidate order.

## Contents

```@contents
Pages = ["reference.md"]
Depth = 3
```

## Core utilities

```@autodocs
Modules = [
    LineCableModels,
    LineCableModels.Grammar,
    LineCableModels.Units,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Data model and materials

```@autodocs
Modules = [
    LineCableModels.Materials,
    LineCableModels.DataModel,
    LineCableModels.DataModel.BaseParams,
    LineCableModels.EarthProps,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Line-parameter engine

```@autodocs
Modules = [
    LineCableModels.Engine,
    LineCableModels.Engine.EarthAdmittance,
    LineCableModels.Engine.EarthProperties,
    LineCableModels.Engine.EarthImpedance,
    LineCableModels.Engine.EHEM,
    LineCableModels.Engine.InsulationAdmittance,
    LineCableModels.Engine.InsulationImpedance,
    LineCableModels.Engine.InternalImpedance,
    LineCableModels.Engine.Transforms,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Parametric and uncertainty modelling

```@autodocs
Modules = [
    LineCableModels.ParametricBuilder,
    LineCableModels.ParametricBuilder.Conductor,
    LineCableModels.ParametricBuilder.Insulator,
    LineCableModels.ParametricBuilder.WirePatterns,
    LineCableModels.UQ,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Plot definitions and recipes

```@autodocs
Modules = [
    LineCableModels.PlotBuilder,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Reports and tables

```@autodocs
Modules = [LineCableModels.ReportBuilder]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Import and export

```@autodocs
Modules = [LineCableModels.ImportExport]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Index

```@index
Pages = ["reference.md"]
Order = [:module, :constant, :type, :function, :macro]
```
