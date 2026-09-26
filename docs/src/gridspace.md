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
ParametricBuilder/UQ  collect stored calculation data
```

`DataModel` constructors validate physical invariants. Engine calculates
core results. UQ performs repeated stochastic realization and aggregation.

## Core invariants

A Gridspace has the semantic shape:

```julia
Gridspace{Target}(build, grids::Tuple; combine=:product)
Gridspace{Target}(grids::Tuple; combine=:product)
```

Every member of `grids` must already be a `Grid` or nested `Gridspace`. The
constructor never interprets a raw tuple, vector, matrix, or other domain value
as an axis. `Target` is the semantic result family used for dispatch and
`build` is the callable that constructs it. A nonempty deterministic space
advertises a concrete iterator element type when Julia can prove that type
without evaluating a point. Otherwise—including uncertainty-bearing and empty
spaces—its iterator uses `Base.EltypeUnknown`; a Gridspace never advertises a
`UnionAll` element type. `combine` is normalized into the concrete Gridspace
type. The space is lazy and has an analytic length. `rand(space)` selects and
realizes one point without collecting the space.

The internal selected point contains only the callable and its selected
arguments. The point is unresolved, unexported, and temporary. The point never enters a
completed calculation result.

## Grid is the variation marker

`Grid` is the only public marker for finite variation:

```julia
Grid((1.0, 2.0, 3.0))
Grid(1.0)                              # one deterministic point
Grid((1.0, 2.0), (1.0, 5.0))          # nominal × relative error [%]
Grid((1.0, 2.0), AbsoluteError(0.1))   # nominal × absolute error
```

The uncertainty-bearing forms yield `UncertainValue(nominal, sigma)`
descriptors. An `UncertainValue` does not depend on an uncertainty package. The descriptor becomes
Measurements values during direct propagation or ordinary scalars during
Monte Carlo realization.

The constructor laws are:

```text
Grid(existing Grid)       = the existing Grid
Grid(existing Gridspace)  = the existing Gridspace
Grid(tuple or array)      = its elements are alternatives
Grid(other value)         = one alternative
```

Grid instances carry no selection identity. Reusing an instance in two source
positions has the same behavior as placing two equal, separately constructed
Grids in those positions.

## Collections are atomic until explicitly varied

Public domain builders know which of their inputs are complete domain values.
An ordinary tuple, vector, or matrix therefore remains atomic:

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
ordering is preserved, so the first source changes fastest:

```julia
space = Gridspace{Tuple}(
    tuple,
    (Grid((1, 2, 3)), Grid((10, 20))),
)

collect(space)
# [(1, 10), (2, 10), (3, 10), (1, 20), (2, 20), (3, 20)]
```

Its length is the product of the direct source lengths. Computing `length`
does not traverse or materialize any point.

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

A nested Gridspace is one finite source at its parent. Its selected value stays
unresolved until recursive materialization or realization reaches it. Nested
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
collection member; the collection itself is reconstructed before the scalar
action runs.

The current behavior is:

| Entry point | Scalar or complete input | Explicit varying input |
|---|---|---|
| `Material` | `Materials.Material` | `Gridspace{Materials.Material}` |
| `Disk`, `Shell`, and other physical geometry | concrete primitive or contextual shell | `Gridspace{Primitive}` or `Gridspace{Shell}` |
| `Region`, `Stack`, `Group`, `Assembly`, `Enclosure` | concrete cable part | `Gridspace{TargetPart}` |
| `CableDesign` | `DataModel.CableDesign` | `Gridspace{DataModel.CableDesign}` |
| `at`, `trefoil`, `hflat`, `vflat` | `Pose2` or `Vector{Pose2}` | corresponding `Gridspace` |
| `homogeneous` | `EarthModel` | `Gridspace{EarthModel}` |
| `LineCableSystem` | `DataModel.LineCableSystem` | `Gridspace{DataModel.LineCableSystem}` |
| `LineParametersProblem` | `Engine.LineParametersProblem` | `Gridspace{Engine.LineParametersProblem}` |
| `CableConstantsProblem` | `Engine.CableConstantsProblem` | `Gridspace{Engine.CableConstantsProblem}` |
| `CableConstants` | `Engine.CableConstants` | `Gridspace{Engine.CableConstants}` |
| `Formulation` | `Engine.LineParametersFormulation` | `Gridspace{Engine.LineParametersFormulation}` |
| `CableConstantsFormulation` | `Engine.CableConstantsFormulation` | `Gridspace{Engine.CableConstantsFormulation}` |
| `ModalAnalysisFormulation` | `ModalAnalysis.ModalAnalysisFormulation` | `Gridspace{ModalAnalysis.ModalAnalysisFormulation}` |
| backend formulation constructor | completed backend formulation | target-bearing formulation `Gridspace` |
| `@gridspace` keyword constructor | strict struct | `Gridspace{Target}` |

Scalar-complete calls invoke their domain action immediately. A varying call
stores only that callable and explicit finite sources; each selected point
invokes the same scalar action.

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

Avoid `Function`-typed fields in frequently called builders. A captured closure is suitable
for local experiments. Reusable callables use concrete types.
Do not introduce a passive record merely to store arguments that another
function immediately unpacks.

`@gridspace` applies the same rule to a keyword-constructed struct. The macro
retains the strict positional constructor, returns the struct immediately for
scalar keyword input, and creates a Gridspace only when a field is an explicit
finite source.

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

The result is a `ParametricResult` in all three cases. A scalar problem forms
a singleton problem axis. `options` belong to the core computations; an
existing `ParametricProblem` retains its stored options. For traversal settings
such as retaining supplemental details, select `Combinatorial` explicitly:

```julia
run = compute(ParametricProblem(problem_space, ComputationOptions(options)),
    Combinatorial(formulation_space; options=(retain_details=true,)))
```

`Combinatorial` accepts one completed formulation, a deterministic
target-bearing formulation `Gridspace`, or a deterministic `Grid` containing
completed formulations. Formulation points are resolved once. Traversal then
materializes each selected problem exactly once and gives the complete
formulation vector to `compute`:

```text
resolve every formulation point once

for each `Gridpoint{Problem}`
    materialise the scalar problem once
    compute(problem, resolved_formulations)
end
```

Every problem/formulation pair is evaluated. Composition inside a formulation
constructor remains local to that formulation: `combine=:product` or
`combine=:zip` determines its formulation points, while the outer
problem/formulation relation is always Cartesian.

The generic vector `compute` method delegates to established scalar dispatch.
Owners may specialize that vector method to share immutable lowering work. The
Coaxial line-parameter path validates and flattens every selected design once
per problem point, constructs its formulation-independent local input once,
then allocates and solves one independent workspace per formulation. The
cable-constant path independently follows the same one-flatten rule. Mutable
matrices, formula-dependent earth data, reduction maps, and diagnostic storage
are never shared. Modal transformation needs no special lowering and uses the
generic route on the same phase-domain matrices.

`ParametricResult` retains both axes:

```julia
run.axes.problems
run.axes.formulations
run[problem_index, formulation_index]
```

Its `values` vector and ordinary linear iteration remain available. Storage is
column-major in `(problem, formulation)` coordinates: the problem index varies
fastest. Thus every selected value can be traced to its completed formulation,
for example with
`formula_id(run.axes.formulations[j].methods.earth_impedance)`. The result does
not retain unresolved points or traversal state.

Direct linear propagation uses the same traversal. The Measurements extension
changes only how an uncertain descriptor materializes. `LinearErrorResult`
likewise stores only its formulation and ordered core results.

ParametricBuilder owns this shared traversal as the qualified `traverse`
method. It computes the first problem/formulation batch, allocates a vector of
that exact result type with the analytic Cartesian cardinality, and rejects any
later type change. Optional detail records follow the same rule and are
resolved through each scalar formulation's computation owner. `Combinatorial`
constructs its result space from `values`, `axes`, and `details`. `LinearError`
continues to use one scalar formulation and consumes `values` and `details`;
Monte Carlo does not use this traversal.

Monte Carlo selects each outer point once, derives a deterministic point seed,
and repeatedly realizes that same point:

```text
for each selected outer point
    until the requested successful trials are collected
        redraw uncertain leaves within that point
        build a fresh complete core problem
        Engine.compute
        optionally reject DomainError realisations under a bounded retry settings
    aggregate that point's draws
end
```

Multiple nominal/error points therefore produce multiple aggregates, not one
mixture. `MonteCarloResult` directly owns sample-mean core results, statistics,
optional retained samples, optional histograms, the root seed, point seeds,
and trial counts.

The default `on_error=:fail` settings rethrows every exception. With
`options=(retain_details=true, on_error=:retry, max_failures=n)`, only
`DomainError` is treated as an unsupported realization. Rejected draws do not
enter samples or statistics, and retry stops when the requested accepted-trial
count is reached or `n` failures have occurred. This estimates the conditional
distribution of the output given that problem construction and computation
succeed; the retained failure summary makes the conditioning rate explicit.

For cable-constant Monte Carlo calculations, the representative stored in the
result space remains a `CableConstants` core result. Retained samples,
statistics, and histograms are concrete named tuples with keys `R`, `L`, `C`,
and `G`; cable samples have assembly × trial dimensions. Line-parameter
products retain conductor × conductor × frequency × trial dimensions. These
are internal point-aligned storage, not result types or observation surfaces.
`MonteCarloResult` validates their keys and dimensions and owns every public
observation method.

### Observe and compare retained uncertainty

Statistical selectors use the same observation grammar as deterministic values:

```julia
using Statistics
using LineCableModels.ReportBuilder: BenchmarkTableDefinition
request = @observe (statistics, R, mean)[1, :, :, :]
table = observables(mc_result, (request, (statistics, R, std, 1)); length_unit=:base)
```

Both MC and LEP expose selected `mean` and `std`. MC means are empirical;
LEP means are first-order nominal predictions. `std` is physical spread, not
uncertainty of the estimated mean. A full MC request `(statistics, R, 1)`
publishes mean/std/min/q05/median/q95/max and the successful trial count.
`Base.Fix2(quantile, 0.05)` selects the retained fifth percentile; unsupported
percentiles fail rather than interpolate an invented distribution.
X/B scale the retained L/C summaries by positive `2πf`. Complex Z/Y support
mean and the nonnegative complex standard deviation; ordered complex
percentiles are undefined. Joint samples and histograms remain separate
products, available only when retained or derivable from retained samples.

```julia
definition = BenchmarkTableDefinition(((statistics, R, mean), (statistics, R, std));
    bands=(:all, :dc, :harmonic, :narrow, :wide))
comparison = report(definition, (reference=mc_result, candidate=lep_result))
comparison.table.features   # Numeric formulation rows, frequency-band columns
comparison.table.statistics
comparison.table.sampling
```

Results without terminal identities require explicit `(result, metadata)`
operands. Multiple reference points require explicit `pairing`; populations
are never pooled. The shared Engine RMS applies the same two-sided numerical
resolution rule to each selected statistic and frequency band.

`confidence(mc_result, point)` reports the actual simultaneous DKW bound,
configured target, marginal count, trial count and conditioning. It counts
both matrix orientations conservatively and applies per outer point, not
campaign-wide. Fixed trials need not certify the configured target. Its
`mean_standard_error` is `s/sqrt(n)` for `n>1`, not an exact confidence interval;
one trial cannot establish population spread. No LEP distribution is inferred
from two moments, and no report initiates sampling or physical computation.

Portable scientific records preserve all retained MC products and shared MC/LEP
Measurement sources, including signed sensitivities. Built-in `:normal` and
`:uniform` law choices are retained directly; Distributions.jl `Normal` and
`Uniform` records retain their actual parameters. Other custom input laws or
formulation declarations require an owner-provided codec and fail explicitly
when no codec exists. They are never changed into a supported law on recovery.

All completed result spaces are one-dimensional finite Julia collections.
Iteration and indexing return one stored core result per calculation. A
`ParametricResult` with several formulations contains the Cartesian
problem/formulation cardinality in its documented storage order. Monte Carlo
iteration returns the stored uncertainty-bearing
core result constructed during aggregation from each point's sample means and
sample standard deviations; individual trials
remain available only through `samples`. Standard `first`, `last`, `only`,
`collect`, `map`, and `zip` operations apply. `only` asserts singleton
cardinality and performs no statistical selection or result transport.

## Transporting completed result spaces

A completed result space enters another scalar calculation through the target
`Gridspace` constructor:

```julia
modal_problems = Gridspace{ModalAnalysisProblem}(phase_results)
```

For example, starting with two completed phase scans and a completed line
formulation, the full modal and finite chain uses public collection constructors
and compute methods throughout:

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
elements; `segments` stores concrete scalar `PropagationParameters` elements.
The first axis of each completed result space follows the source order; for
the product, the linear index is `source_index + (formulation_index-1)*N`.
The quantity arrays
remain single outer elements under ordinary Julia broadcast. A completed
formulation `Grid` can be used directly; it is never wrapped as a second
modal formulation.

The dispatched method `Gridspace{Target}(source::SourceResult)` defines how the
source result family supplies arguments to the target problem constructor. The
transport preserves source cardinality and order unless that source-specific
method explicitly documents another operation. It produces a target-bearing
problem space, never a nested result envelope.

`ParametricResult` transports its completed combinatorial results directly.
`LinearErrorResult` and `MonteCarloResult` transport their stored
uncertainty-bearing results directly. Load `Measurements` before MC computation;
the extension is required before sampling starts. Aggregation constructs each
marginal from the accepted raw samples, independently of histogram bins and
sample-retention options. The uncertainty is the output sample standard
deviation, not the standard error of its mean. This marginal surrogate does not
recover covariance between observables. Repeated indexing, `uncertain`, and
transport reuse its existing uncertainty-source identities without rebuilding
Measurements.

Only result families with defined semantics are admitted. An unsupported
source/target pair raises an error naming both types and directs the caller to
`?Gridspace` and this page. Extension code adds a transport by defining:

```julia
import LineCableModels: Gridspace

function Gridspace{Target}(source::OwnedResultSpace)
    # Expose stored source values, then return Gridspace{Target}.
end
```

This keeps every scientific calculation scalar. Finite composition recurses
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
positions. By contrast, two separately declared uncertain source positions are
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
clearance is positive; shared positive scaling preserves that sign and the
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
Normal sampling remains the default and has unbounded support; finite reserve
distances cannot make all independent normal draws feasible. Retry mode
estimates a distribution conditional on success, not the original input law.

Executable native checkpoints retain the joint builder and its sources.
Marginal Measurement JSON retains values and standard uncertainties only; it is
not a format for archiving the joint statistical law or covariance.

## Optional package extensions

The core package declares uncertainty without loading Measurements or
Distributions.

- Loading Measurements adds direct materialization of `UncertainValue` while
  retaining exact structural reuse.
- Loading Distributions adds standardised univariate sampling families. The
  selected distribution must have finite mean and positive finite standard
  deviation. Samples are transformed to the descriptor's nominal value and
  standard uncertainty.

Neither extension knows about finite-source identity, result presentation, or
Engine internals.

## Performance and conformance

The implementation relies on tuple-specialized recursion and Julia's public
product and zip iterators. The implementation guarantees:

- `length` is analytic and never enumerates points.
- Product and zip traversal are linear in yielded work.
- Materialization and realization use no dictionary or identity lookup.
- Immutable scalar targets infer through selection, materialization, and
  realization.
- After warmup, deterministic iteration and a bare 10,000-realization scalar
  loop add zero heap allocations.
- Full cable and line construction may allocate only what their existing
  vectors, domain constructors, Engine computation, and requested result
  storage intrinsically require.

The conformance suite in `test/unit/parametricbuilder/conformance.jl` checks
these properties, exact structural reuse, explicit variation, scalar public
construction, and the absence of a random-access Gridspace API.

## Implementation map

The implementation is split across:

- `src/grid.jl`: finite values and uncertainty descriptors.
- `src/gridspace.jl`: composition, point selection, recursive
  materialization, and realization.
- `src/parametricbuilder/macros.jl`: strict scalar construction and explicit
  Gridspace lifting for `@gridspace`.
- material, cable, position, and system files: scalar construction and lifting
  rules and concrete callable algorithms.
- `src/parametricbuilder/traversal.jl`: combinatorial traversal.
- `src/modalanalysis/delegation.jl`: public phase-to-modal composition;
  `src/modalanalysis/problems.jl`, `compute.jl`, and `propagation.jl` own
  modal scalar computation and finite segment binding.
- `src/uq/linearerror.jl` and `src/uq/montecarlo/compute.jl`: direct and
  repeated stochastic traversal.
- Measurements and Distributions extensions: dependency-specific uncertainty
  behavior only.
