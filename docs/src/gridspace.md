# Gridspace

Gridspace combines finite sets of inputs and constructs one object for each
selected combination. The cable API uses it to construct designs, systems,
and problems with varying parameters.

The calculation sequence is:

```text
Grid                 declares explicit finite variation
Gridspace            composes finite sources and selects one point
callable              invokes the existing complete construction action
Engine.compute        evaluates one complete core problem
ParametricBuilder/UQ  collect stored computation data
```

`DataModel` constructors validate physical invariants. Engine calculates
core results. UQ performs repeated stochastic realization and aggregation.

## Core invariants

A Gridspace has this structure:

```julia
Gridspace{Target}(build, grids::Tuple; combine=:product)
Gridspace{Target}(grids::Tuple; combine=:product)
```

Every member of `grids` must already be a `Grid` or nested `Gridspace`. Wrap
finite alternatives in `Grid` before passing them to this constructor.
`Target` is the result family used for dispatch and
`build` is the callable that constructs it. `Target` is the element type of the space,
and its iterator uses `Base.EltypeUnknown`: each point has the concrete type that `build`
returns. `combine` selects `Product` or `Zip`, which is part of the Gridspace type. A
zip pairs the direct parameters of one call, so it refuses a nested space that varies. A
nested space of one point pairs as that point. To pair complete objects, pass a `Grid`
of them. The space is lazy and has an analytic length. `rand(space)` selects and
realizes one point without collecting the space.

## Grid is the variation marker

`Grid` is the only public marker for finite variation:

```julia
Grid((1.0, 2.0, 3.0))
Grid(1.0)                              # one deterministic point
Grid((1.0, 2.0), (1.0, 5.0))          # nominal × relative error [%]
Grid((1.0, 2.0), AbsoluteError(0.1))   # nominal × absolute error
```

The uncertainty-bearing forms yield `UncertainValue(nominal, sigma)`
descriptors defined by the core package. A descriptor becomes a Measurements
value during direct propagation or an ordinary scalar during Monte Carlo realization.

The constructor laws are:

```text
Grid(existing Grid)       = the existing Grid
Grid(existing Gridspace)  = the existing Gridspace
Grid(tuple or array)      = its elements are alternatives
Grid(other value)         = one alternative
```

Grid instances include no selection identity. Reusing an instance in two source
positions has the same behavior as placing two equal, separately constructed
Grids in those positions.

## Collections are atomic until explicitly varied

Public domain builders know which of their inputs are complete domain values.
An ordinary tuple, vector, or matrix remains atomic:

```julia
frequencies = [50.0, 100.0, 1000.0]       # one frequency scan
payload = [1.0 2.0; 3.0 4.0]              # one matrix

frequency_sets = Grid((
    [50.0, 100.0],
    [50.0, 500.0, 5_000.0],
))                                         # two complete scans
```

When one sibling varies, the builder wraps every ordinary sibling in a
one-point source.
The matrix or frequency vector remains one complete value at every point.
Calling `Grid(matrix)` is different and explicit: the matrix elements become
alternatives because the caller requested that variation.

A tuple or vector containing a `Grid` or `Gridspace` is not atomic: the
explicit finite sources are composed and the collection is rebuilt from their
selected, completed values. This lets domain collections such as
`[design_space, design_space]` reach their scalar builder as completed
`CableDesign` objects. Ordinary collections containing no explicit finite
source remain one complete value.

## Product and zip composition

Composition is local to each Gridspace node.

### Product

`:product` is the default and forms the lazy Cartesian product. Julia product
ordering is preserved, with the first source changing fastest:

```julia
space = Gridspace{Tuple}(
    tuple,
    (Grid((1, 2, 3)), Grid((10, 20))),
)

collect(space)
# [(1, 10), (2, 10), (3, 10), (1, 20), (2, 20), (3, 20)]
```

Its length is computed directly as the product of the source lengths.

### Zip

`:zip` pairs non-singleton sources row by row and broadcasts singleton
sources:

```julia
space = Gridspace{Tuple}(
    tuple,
    (Grid((1, 2, 3)), Grid((10, 20, 30)), Grid(:fixed));
    combine=:zip,
)

collect(space)
# [(1, 10, :fixed), (2, 20, :fixed), (3, 30, :fixed)]
```

All non-singleton direct sources must have equal cardinality. A mismatch
throws `DimensionMismatch` when the Gridspace is constructed, before any point
is materialized. Zip traversal is linear in the number of rows.

### Nesting

A nested Gridspace is one finite source at its parent. Its selected value remains unresolved until recursive materialization or realization reaches it. Nested
resolution lets one child zip local parameters while its parent forms a
Cartesian product:

```julia
paired = Gridspace{Tuple}(
    tuple,
    (Grid((1, 2)), Grid((10, 20)));
    combine=:zip,
)

outer = Gridspace{Tuple}(tuple, (paired, Grid((:a, :b))))
collect(outer)
# [((1, 10), :a), ((2, 20), :a),
#  ((1, 10), :b), ((2, 20), :b)]
```

## Public constructors

Gridspace delays finite selection. Each lifted public action applies one rule:

```text
no admitted Grid or Gridspace       -> invoke the scalar action now
at least one explicit source        -> preserve sources, singleton-wrap siblings,
                                       and return a Gridspace
```

For domain arguments that are tuples or vectors, an admitted source may be a
collection member. The collection itself is reconstructed before the scalar
action runs.

The current behavior is:

| Entry point | Scalar or complete input | Explicit varying input |
|---|---|---|
| `Material` | `Materials.Material` | `Gridspace{Materials.Material}` |
| `Disk`, `Shell`, and other physical geometry | concrete primitive or contextual shell | `Gridspace{Primitive}` or `Gridspace{Shell}` |
| `Region`, `Stack`, `Group`, `Assembly`, `Enclosure` | concrete cable part | `Gridspace{TargetPart}` |
| `CableDesign` | `DataModel.CableDesign` | `Gridspace{DataModel.CableDesign}` |
| `at`, `trefoil`, `hflat`, `vflat` | `Pose2`, placed member, placed-design record or a formation | corresponding `Gridspace` |
| `layers`, `terminal`, `solid`, `core`, `insulation` and the other builders | the declared part | `Gridspace` of the builder's target |
| `formula` | `FormulaDefinition` | `Gridspace{FormulaDefinition}` |
| `homogeneous` | `EarthModel` | `Gridspace{EarthModel}` |
| `LineCableSystem` | `DataModel.LineCableSystem` | `Gridspace{DataModel.LineCableSystem}` |
| `LineParametersProblem` | `Engine.LineParametersProblem` | `Gridspace{Engine.LineParametersProblem}` |
| `CableConstantsProblem` | `Engine.CableConstantsProblem` | `Gridspace{Engine.CableConstantsProblem}` |
| `CableConstants` | `Engine.CableConstants` | `ParametricResult` from `compute` |
| `Formulation` | `Engine.LineParametersFormulation` | `Gridspace{Engine.LineParametersFormulation}` |
| `CableConstantsFormulation` | `Engine.CableConstantsFormulation` | `Gridspace{Engine.CableConstantsFormulation}` |
| `ModalAnalysisFormulation` | `ModalAnalysis.ModalAnalysisFormulation` | `Gridspace{ModalAnalysis.ModalAnalysisFormulation}` |
| backend formulation constructor | completed backend formulation | target-bearing formulation `Gridspace` |
| hand-written lift of a user type | the user's type | `Gridspace` of the user's type |

Scalar-complete calls invoke their domain action immediately. A varying call
stores only that callable and explicit finite sources. Each selected point
invokes the same scalar action.

Each builder lifts over all its arguments and keywords at its own level, and its points build
their inner objects from scalars. In a product, the first source varies fastest, in the order
the call writes them. `wires(cu; shape, n = 6, lay = Grid(...), r = Grid(...))` varies `lay`
fastest. The [developers guide](developers.md#Lifting) shows how to define a builder and how
to extend one for a type of your own.

## Gridspace callables

A Gridspace callable should be a concrete immutable functor whose field types
are concrete. The callable should perform one construction step and delegate physical
validation to the target constructors:

```julia
struct PairValue end
(::PairValue)(left, right) = (left, right)

space = Gridspace{Tuple}(
    PairValue(),
    (Grid((1, 2)), Grid((10, 20))),
)
```

Avoid `Function`-typed fields in frequently called builders. A function that captures
local variables is suitable for local experiments. Reusable callables use concrete types.
Do not introduce a passive record merely to store arguments that another
function immediately unpacks.

A user type lifts with one method that calls `parameterize`, as the package's constructors
do. Each position accepts its scalar type or a `Grid` or `Gridspace`, and the type is the
caller:

```julia
struct LabeledValue{T <: Real}
    value::T
    label::Symbol
end

function LabeledValue(
        value::Union{Real, AbstractGrid, Gridspace},
        label::Union{Symbol, AbstractGrid, Gridspace} = :default;
        combine::Symbol = :product
)
    return LineCableModels.parameterize(
        LabeledValue, LabeledValue, (value, label); combine)
end

LabeledValue(1.0)                    # LabeledValue(1.0, :default)
LabeledValue(Grid((1.0, 2.0)), :ok)  # a Gridspace with two points
```

The constructor `LabeledValue(value::T, label::Symbol)` is more specific than the lift, so
`parameterize` reaches it for scalar inputs and at each point of the space.

## Materialization and realization

Ordinary iteration selects an internal target-bearing `Gridpoint{Target}` and recursively
materializes its arguments. Deterministic values pass through unchanged.
Nested points invoke their own callable before the parent callable is invoked.
After loading Measurements, an `UncertainValue` materializes as one
`Measurement`.

Stochastic realization follows the same recursion with a caller-owned random
number generator. Only `UncertainValue` leaves are redrawn. Deterministic
selections remain fixed. A zero-sigma descriptor resolves deterministically.

Materialization and realization are internal. Developers extend the public grammar through
concrete builders and supported uncertainty extensions, not
by exposing unresolved points as application data.

## Higher-order computation

Passing a formulation `Gridspace` directly to `compute` selects default
`Combinatorial` traversal:

```julia
run = compute(problem, formulation_space; options=(;))
run = compute(problem_space, formulation_space; options=(;))
run = compute(ParametricProblem(problem_space, ComputationOptions(options)), formulation_space)
```

The result is a `ParametricResult` in all three cases: the first two calls compute the
`ParametricProblem` of their problem and `options`, as the third does. A scalar problem
forms a singleton problem axis. `options` belong to the core computations. An
existing `ParametricProblem` retains its stored options. For traversal settings
such as retaining supplemental details, select `Combinatorial` explicitly:

```julia
run = compute(ParametricProblem(problem_space, ComputationOptions(options)),
    Combinatorial(formulation_space; options=(retain_details=true,)))
```

`Combinatorial` accepts one completed formulation, a deterministic
target-bearing formulation `Gridspace`, or a deterministic `Grid` containing
completed formulations. A formulation space contains the formulations of one backend,
whose formula choices may vary. Comparing backends takes one computation per backend. The runner, `traverse`, materializes and validates the
problems and formulations before the first point, then computes problem by problem:

```text
materialize and validate every problem and formulation

for each problem
    resolve(problem, formulations; options)      # once per problem
    for each formulation
        compute(resolved, formulation; options)  # one point
    end
end
```

Each problem-formulation pair is evaluated. Composition inside a formulation
constructor remains local to that formulation: `combine=:product` or
`combine=:zip` determines its formulation points, while the outer
problem-formulation relation is always Cartesian.

`resolve` shares the work that the formulations of one problem need equally. The
coaxial line-parameter path flattens the selected designs and builds the
formulation-independent line input, both once per problem. It completes the internal
potential coefficients of the formulations together. Each point then allocates and
solves one independent workspace. The cable-constant path follows the same
one-flatten rule. Mutable matrices, formula-dependent earth data, reduction maps,
and diagnostic storage belong to each workspace. Modal analysis uses the default `resolve`.
PSCAD and FEM share one table of native runs per problem, so a later point with the
same native input reuses the run that an earlier point launched.

`ParametricResult` retains both axes:

```julia
run.axes.problems
run.axes.formulations
run[problem_index, formulation_index]
```

Its `values` vector and ordinary linear iteration remain available. Storage is
column-major in `(problem, formulation)` coordinates: the problem index varies
fastest. The formulation for column `j` is stored in `run.axes.formulations[j]`.
Use `formula_id(run.axes.formulations[j].formulas.earth_impedance)` to identify
its earth-impedance method.

Direct linear propagation uses the same traversal. The Measurements extension
changes only how an uncertain descriptor materializes. `LinearErrorResult`
stores only its formulation and ordered core results.

`traverse` computes the first point and allocates a vector of that exact result type
with the analytic Cartesian cardinality. It writes every later value at its index in
the space. A scalar problem and a scalar formulation are a run of one point, and
`compute` returns its value. Optional detail records are resolved through each
formulation's computation owner. `Combinatorial` keeps its result space as `values`,
`axes` and `details`. `LinearError` uses one scalar formulation and consumes
`values` and `details`.

Monte Carlo selects each outer point once and derives a deterministic seed for that point. It repeatedly realizes the same point:

```text
for each selected outer point
    until the requested successful trials are collected
        redraw uncertain leaves within that point
        build a fresh complete core problem
        Engine.compute
        optionally reject DomainError realizations with a bounded number of resampling attempts
    aggregate that point's draws
end
```

Each nominal and error point produces its own aggregate. `MonteCarloResult` directly owns sample-mean core results, statistics,
optional retained samples, optional histograms, the root seed, point seeds,
and trial counts.

The default `on_error=:fail` rethrows every exception. With
`options=(retain_details=true, on_error=:resample, max_failures=n)`, only
`DomainError` is treated as an unsupported realization. Rejected draws do not
enter samples or statistics, and resampling stops when the requested accepted-trial
count is reached or `n` failures have occurred. This estimates the conditional
distribution of the output when problem construction and computation
succeed. The retained failure summary makes the conditioning rate explicit.

For cable-constant Monte Carlo computations, the representative stored in the
result space remains a `CableConstants` core result. Retained samples,
statistics, and histograms are concrete named tuples with keys `R`, `L`, `C`,
and `G`. Cable samples have assembly × trial dimensions. Line-parameter
products retain conductor × conductor × frequency × trial dimensions.
`MonteCarloResult` validates these point-aligned arrays' keys and dimensions
and supplies their public observation methods.

### Observe and compare retained uncertainty

Statistical selectors use the same observation grammar as deterministic values:

```julia
using Statistics
using LineCableModels.ReportBuilder: BenchmarkTableDefinition
request = @observe (statistics, R, mean)[1, :, :, :]
table = observables(mc_result, (request, (statistics, R, std, 1)); length_unit=:base)
```

Both MC and LEP expose selected `mean` and `std`. MC means are empirical.
LEP means are first-order nominal predictions. `std` is physical spread, not
uncertainty of the estimated mean. A full MC request `(statistics, R, 1)`
publishes mean, standard deviation, minimum, fifth percentile, median, 95th percentile, maximum and the successful trial count.
`Base.Fix2(quantile, 0.05)` selects the retained fifth percentile. Unsupported
percentiles fail rather than interpolate an invented distribution.
X/B scale the retained L/C summaries by positive `2πf`. Complex Z/Y support
mean and the nonnegative complex standard deviation. Ordered complex
percentiles are undefined. Joint samples and histograms remain separate
products, available only when retained or derivable from retained samples.

```julia
definition = BenchmarkTableDefinition(((statistics, R, mean), (statistics, R, std));
    bands=(:all, :dc, :harmonic, :narrow, :wide))
comparison = report(definition, (reference=mc_result, result=lep_result))
comparison.table.features   # Numeric formulation rows, frequency-band columns
comparison.table.statistics
comparison.table.sampling
```

Results without terminal identities require explicit `(result, metadata)`
operands. Multiple reference points require explicit `pairing`. Populations
are never pooled. The shared Engine RMS applies the same two-sided numerical
resolution rule to each selected statistic and frequency band.

`confidence(mc_result, point)` reports the actual simultaneous DKW bound,
configured target, marginal count, trial count and conditioning. It counts
both matrix orientations conservatively and applies per outer point, not
campaign-wide. Fixed trials need not certify the configured target. Its
`mean_standard_error` is `s/sqrt(n)` for `n>1`, not an exact confidence interval.
One trial cannot establish population spread. No LEP distribution is inferred
from two moments, and no report initiates sampling or physical computation.

A Monte Carlo or linear-error result saved with `save` keeps its retained products, its
formulation, the law it sampled with and its seeds, and `import_data(path, :result)` reads
it with the same code (see [Saving results](@ref)). Values that shared a Measurement source
in the saved result share one again after loading.

All completed result spaces are one-dimensional finite Julia collections.
Iteration and indexing return one stored core result per computation. A
`ParametricResult` with several formulations contains the Cartesian
problem-formulation cardinality in its documented storage order. Monte Carlo
iteration returns the stored uncertainty-bearing
core result constructed during aggregation from each point's sample means and
sample standard deviations. Individual trials
remain available only through `samples`. Standard `first`, `last`, `only`,
`collect`, `map`, and `zip` operations apply. `only` asserts singleton
cardinality and does not perform statistical selection or result transport.

## Transporting completed result spaces

A completed result space enters another scalar computation through the target
`Gridspace` constructor:

```julia
modal_problems = Gridspace{ModalAnalysisProblem}(phase_results)
```

The following example starts with two completed phase scans and a completed
line formulation. It uses public collection constructors and compute methods
to calculate modal parameters and then finite-length responses:

```julia
# phase_a and phase_b are completed phase LineParameters with declared lengths.
phase_results = ParametricResult(Combinatorial(line_formulation),
    [phase_a, phase_b])
modal_formulation = ModalAnalysisFormulation(:default)

# One scalar action and direct target-bearing Gridspace action.
scalar_modal = compute(ModalAnalysisProblem(phase_a), modal_formulation)
modal_problems = Gridspace{ModalAnalysisProblem}(phase_results)
modal_results = compute(modal_problems, modal_formulation)
length(modal_results) == 2                 # N → N, original source order

# Explicit Combinatorial uses the same finite traversal.
explicit_results = compute(ParametricProblem(modal_problems),
    Combinatorial(modal_formulation))
length(explicit_results) == 2

two_formulations = ModalAnalysisFormulation(
    Grid((modal_formulation, modal_formulation)))
product_results = compute(modal_problems, two_formulations)
length(product_results) == 4               # N × M, source index varies fastest

segments = collect(Gridspace{PropagationParameters}(modal_results))
short_segments = PropagationParameters(scalar_modal;
    line_length=Grid((150.0, 300.0)))
roots = gamma.(modal_results)              # each element is a mode × frequency array
voltage_bases = Tv.(modal_results)          # each element is a phase × mode × frequency array
responses = H.(segments)                   # each element is a mode × frequency array
```

`modal_results` and `product_results` store concrete scalar `LineParameters`
elements. `segments` stores concrete scalar `PropagationParameters` elements.
The first axis of each completed result space follows the source order. For
the product, the linear index is `source_index + (formulation_index-1)*N`.
The quantity arrays
remain single outer elements under ordinary Julia broadcast. A completed
formulation `Grid` can be passed directly to `compute`.

The dispatched method `Gridspace{Target}(source::SourceResult)` defines how arguments from the source family are supplied to the target problem constructor. The
transport preserves source cardinality and order unless that source-specific
method explicitly documents another operation. It produces a `Gridspace` of
the target problem type.

`ParametricResult` transports its completed combinatorial results directly.
`LinearErrorResult` and `MonteCarloResult` transport their stored
uncertainty-bearing results directly. Load `Measurements` before MC computation.
The extension is required before sampling starts. Aggregation constructs each
marginal from the accepted raw samples, independently of histogram bins and
sample-retention options. The uncertainty is the output sample standard
deviation, not the standard error of its mean. Covariance between observables remains unavailable from this marginal surrogate. Repeated indexing, `uncertain`, and
transport reuse its existing uncertainty-source identities without rebuilding
Measurements.

Only result families with defined semantics are admitted. An unsupported
source-target pair raises an error identifying both types and directs the caller to
`?Gridspace` and this page. Extension code adds a transport by defining:

```julia
import LineCableModels: Gridspace

function Gridspace{Target}(source::OwnedResultSpace)
    # Expose stored source values, then return Gridspace{Target}.
end
```

This keeps every scientific computation scalar. Finite composition recurses
through the same typed result-to-problem conversion:

```text
ResultSpace{A} → Gridspace{ProblemB} → ResultSpace{B}
```

## Pairing, exact reuse, and correlation

Zip pairing, exact argument reuse, and stochastic correlation have different
semantics.

`combine=:zip` is deterministic row pairing between finite sources. Zip pairing says
nothing about covariance.

Exact reuse is structural. Pass one selected uncertain argument once to a
callable and use that argument more than once:

```julia
struct Duplicate end
(::Duplicate)(value) = (value, value)

space = Gridspace{Tuple}(
    Duplicate(),
    (Grid(10.0, AbsoluteError(0.5)),),
)
```

Direct propagation constructs one Measurement and both tuple positions retain
that variable. Monte Carlo draws once and passes the same scalar to both
positions. But two separately declared uncertain source positions are
independent, even when they contain the same Grid instance or numerically equal
descriptors.

General correlation between distinct variables is a UQ concern. Correlation requires a
joint stochastic source or distribution that returns a tuple or vector sample
consumed by one builder. Gridspace does not infer or register correlation.

### Feasible geometric dependence

For a fixed-count wire ring, varying the ring radius and wire diameter
independently can produce overlaps even when the nominal ring fits. Declare the
intended dependence in an ordinary builder, before constructing physical parts:

```julia
using LineCableModels, Measurements, Random

joint = Gridspace{NamedTuple{(:inner, :diameter, :x)}}(
    (scale, x) -> (inner=0.0679scale, diameter=0.006scale, x=x),
    (Grid(1.0, 10.0), Grid(0.0, AbsoluteError(0.002))),
)
geometry = Gridspace{Tuple}(
    p -> (
        Ring(68; r=p.inner+p.diameter/2),
        Disk(p.diameter/2),
        Annulus(p.inner, p.inner+p.diameter, Pose2(p.x, 0.0)),
    ),
    (joint,),
)
propagated = only(geometry)
sampled = rand(Xoshiro(42), geometry; distribution=:uniform)
```

Lengths are in meters. Here the common scale has mean 1, standard deviation
0.1 and uniform support `[1-sqrt(3)*0.1, 1+sqrt(3)*0.1]`. The nominal chord
clearance is positive. Shared positive scaling preserves that sign and the
68-wire inventory throughout this support. The annulus thickness is derived
from the same diameter, and its independent coordinate retains its own
uncertainty. For a complete stack, derive successive layer radii from positive
thicknesses inside the same builder. Nest the resulting source once in the
`LineParametersProblem` builder, then use it with `ParametricProblem` and
either `LinearError` or `MonteCarlo`.

This is a specified correlated model, not independent manufacturing tolerances.
It preserves each length's nominal mean and 10% standard deviation, but derived
areas scale quadratically: their mean is `1.01` times nominal area. Linear
propagation is local and does not include that second-order mean shift.
Normal sampling remains the default and has unbounded support. Finite reserve
distances cannot make all independent normal draws feasible. Resampling
estimates a distribution conditional on success, not the original input law.

A result saved with `save` keeps the joint builder of its space and the sources of its
values. Marginal Measurement JSON retains values and standard uncertainties only. It is
not a format for archiving the joint statistical law or covariance.

## Optional package extensions

The core package declares uncertainty without loading Measurements or
Distributions.

- Loading Measurements adds direct materialization of `UncertainValue` while
  retaining exact structural reuse.
- Loading Distributions adds standardized univariate sampling families. The
  selected distribution must have finite mean and positive finite standard
  deviation. Samples are transformed to the descriptor's nominal value and
  standard uncertainty.

## Performance and conformance

The implementation relies on tuple-specialized recursion and Julia's public
product and zip iterators. The implementation guarantees:

- `length` is computed analytically from source cardinalities.
- Product and zip traversal are linear in yielded work.
- Materialization and realization do not use a dictionary or identity lookup.
- Immutable scalar targets infer through selection, materialization, and
  realization.
- After warmup, deterministic iteration and a bare 10,000-realization scalar
  loop run without heap allocation.
- Allocations during full cable and line construction are limited to the requirements of existing vectors and domain constructors, the Engine computation and the requested result storage.

The conformance tests in `test/unit/commons/gridspace.jl` and
`test/unit/parametricbuilder/conformance.jl` check these properties, exact structural
reuse, explicit variation, and scalar public construction.

## Implementation map

The implementation is split across:

- `src/commons/grid.jl`: finite values and uncertainty descriptors.
- `src/commons/gridspace.jl`: composition, point selection, recursive
  materialization, and realization.
- material, cable, position, and system files: scalar construction and lifting
  rules and concrete callable algorithms.
- `src/parametricbuilder/traversal.jl`: combinatorial traversal and the
  phase-to-modal composition of a `ParametricProblem`.
- `src/engine/modalanalysis/composition.jl`: public phase-to-modal composition.
  `src/engine/modalanalysis/problems.jl`, `compute.jl`, and `propagation.jl` own
  modal scalar computation and line segment binding.
- `src/uq/linearerror.jl` and `src/uq/montecarlo/compute.jl`: direct and
  repeated stochastic traversal.
- Measurements and Distributions extensions: dependency-specific uncertainty
  behavior only.
