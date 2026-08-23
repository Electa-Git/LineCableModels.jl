# API reference

## Modeling and execution grammar

The top-level namespace exposes the declarative modeling API. `Grid` is the
only parameter-variation syntax, `Gridspace{Target}` is a lazy space of complete
targets, and `Gridspace(definition)` is the explicit materialization boundary
for cable, position, earth, and system definitions.

Ordinary collections remain ordinary constructor values. Wrap a collection in
`Grid` only when it is intended to vary. Reusing the same `Grid` deliberately
couples the selection and stochastic realization wherever that grid appears.

```julia
earth = Earth(
    rho=Grid((10.0, 100.0, 1000.0)),
    eps_r=Grid((100.0, 10.0, 5.0));
    combine=:zip,
)
```

Composition is local to each Gridspace node. The example above admits exactly
three earth configurations after `Gridspace(earth)` is constructed; the earth
definition can still participate in a Cartesian product at a parent problem
node.

All numerical work follows `problem → formulation → compute`:

```julia
compute(problem, Formulation())

space = Gridspace(problem_definition)
compute(ParametricProblem(space), Combinatorial(Formulation()))
compute(ParametricProblem(space), LinearError(Formulation()))
compute(
    ParametricProblem(space),
    MonteCarlo(Formulation(); trials=1000, cdf_tol=0.02),
)
```

The higher-order formulation is always explicit. `Combinatorial` evaluates
every admitted configuration, `LinearError` selects established direct
propagation, and `MonteCarlo` samples realizations within each selected outer
configuration. Every cardinality, including one, returns the corresponding
composite result family.

`Formulation` owns physics and numerical-method choices. The separate
computation-options named tuple is shared unchanged by ordinary,
combinatorial, and uncertainty calculations. A materialized line system can
request total rather than per-length Z/Y matrices without altering its
formulation:

```julia
compute(problem, formulation; options=(output_basis=:total,))
```

`output_basis=:total` scales both impedance and admittance by the materialized
system length. `CableConstantsProblem` has no line length and therefore accepts
only the default `:per_length` basis.

## Results and provenance

`CableConstants` stores R/L/C values per metre. `LineParameters` stores
frequency-dependent Z/Y matrices with their domain and `:per_length` or
`:total` basis.

`observables(result)` publishes the stable scientific quantities used by
tables and plots as an immutable named tuple. Array-valued entries use detached
storage, so presentation code cannot mutate the completed result. The primary
keys are:

- `CableConstants`: `:resistance`, `:inductance`, and `:capacitance`;
- `LineParameters`: `:frequency`, `:series_impedance`, and
  `:shunt_admittance`;
- `LineParametersBenchmark`: full series/shunt absolute and relative RMS-error
  keys.

The mathematical `Z`, `Y`, `R`, `X`, `L`, `G`, `B`, and `C` accessors remain
available for explicit derived selection. Plot-only component spellings do not
belong to the common observable grammar.

A complete deterministic traversal returns `ParametricResult{T}`. Direct
propagation returns `LinearErrorResult{T}`, and conditional sampling returns
`MonteCarloResult{T}`. In every family, `T` is the primitive `CableConstants`
or `LineParameters` result rather than another composite result.

Use `result`, `statistics`, `samples`, `histograms`, `uncertain_value`, and
`manifest` to inspect results. Every result details dictionary has the keys
`:failures`, `:samples`, `:histograms`, `:random`, and `:manifest`. A
calculation manifest contains a stable SHA-256 hash over the resolved
parameterization, original problem assumptions, formulation, solver identity,
execution policy, and calculation options.

`DataFrame(result::MonteCarloResult)` renders stored marginal summaries
without repeating the calculation. Cable-constant results produce one R/L/C
table; line-parameter results produce one R/L/C/G table for every matrix entry
and frequency. The displayed `confidence` and `cdf_tol` values describe the DKW
bound below and are not confidence intervals for the sample mean.

After loading a Makie package, retained samples and histograms can be displayed
through the maintained Monte Carlo recipe:

```julia
using CairoMakie

plot(result, :R; mode=:hist, data=:both)
plot(result, :R; mode=:pdf)
plot(result, :R; mode=:ecdf, data=:both)
plot(result, :R; mode=:qq)

# Select one line-parameter matrix entry and frequency:
plot(line_result, :L; ijk=(1, 1, 3), mode=:hist)
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

## Optional uncertainty integrations

The core Grid/Gridspace grammar does not load Measurements.jl or
Distributions.jl. Loading Measurements enables direct `LinearError(inner)`
propagation while preserving shared variables and correlations.
`MonteCarloResult` has no implicit conversion to Measurements values.
Empirical reconstruction is an explicit consumer operation and requires
retained joint samples. Loading Distributions enables compatible univariate
distributions as Monte Carlo samplers and `pdf`/`cdf` evaluation of
`HistogramDensity`. A supplied distribution is standardized to the Grid
descriptor's nominal value and standard deviation, so it must have a finite
mean and a positive finite standard deviation.

## Strict materialized API

The eager cable and system model remains the materialized boundary but is not
exported at top level. Import those types explicitly when strict construction
is required. `Material`, `MaterialsLibrary`, and `CablesLibrary` are public
because they are shared by both modeling modes:

```julia
using LineCableModels: Material
import LineCableModels.DataModel
import LineCableModels.EarthProps

copper = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
core = DataModel.Tubular(0.0, 10e-3, copper)
```

Materialized radial constructors accept resolved numeric inner and outer radii.
Radius-versus-thickness intent and repetition belong to builder definitions.
Group `add!` accepts a fully materialized part of the same scalar type as its
destination.

All eager cable quantities are evaluated at the common material reference
temperature. The materialized `CableDesign` therefore rejects layers with
different `T0` values. Operating temperature belongs only to
`LineParametersProblem`; `SystemBuilder` exposes that setting and [`compute`](@ref)
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
A physically valid search that cannot meet every requested limit still returns ranked
best-effort candidates and concise reasons:

```julia
estimate = make_stranded(1000.0)
closest = estimate[Val(:match)]
fewest_layers = estimate[Val(:layers)]
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
    LineCableModels.Computation,
    LineCableModels.UnitHandler,
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

## Parametric modeling

```@autodocs
Modules = [
    LineCableModels.ParametricBuilder,
    LineCableModels.ParametricBuilder.WirePatterns,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Plot definitions and recipes

```@autodocs
Modules = [
    LineCableModels.PlotBuilder,
    LineCableModels.PlotBuilder.BackendHandler,
]
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
