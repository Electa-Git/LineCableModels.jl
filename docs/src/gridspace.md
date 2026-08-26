# Gridspace

Gridspace combines explicit finite sources, selects one unresolved point, and
invokes a typed callable after resolving that point. The declarative cable API
uses it to delay construction until each finite selection is known.
Gridspace does not compute line parameters, aggregate uncertainty, validate physical
objects independently of their constructors, or store traversal state in
completed results.

The calculation sequence is:

```text
Grid                 declares explicit finite variation
Gridspace            composes finite sources and selects one point
callable builder      constructs the existing eager domain object
Engine.compute        evaluates one complete core problem
ParametricBuilder/UQ  collect stored calculation data
```

`DataModel` constructors validate physical invariants. Engine calculates
core results. UQ performs repeated stochastic realisation and aggregation.

## Core invariants

A Gridspace has the semantic shape:

```julia
Gridspace{Target}(build, grids::Tuple; combine=:product)
Gridspace{Target}(grids::Tuple; combine=:product)
```

Every member of `grids` must already be a `Grid` or nested `Gridspace`. The
constructor never interprets a raw tuple, vector, matrix, or other domain value
as an axis. `Target` is the type promised by iteration, `build` is the callable
that constructs it, and `combine` is normalised into the concrete Gridspace
type. The space is lazy, has an analytic length, and has no random-access or
random-sampling shortcut.

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
Monte Carlo realisation.

The constructor laws are:

```text
Grid(existing Grid)       = the existing Grid
Grid(existing Gridspace)  = the existing Gridspace
Grid(tuple or array)      = its elements are alternatives
Grid(other value)         = one alternative
```

Grid instances carry no selection identity. Reusing an instance in two source
positions has the same behaviour as placing two equal, separately constructed
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
does not traverse or materialise any point.

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
is materialised. Zip traversal is linear in the number of rows.

### Nesting

A nested Gridspace is one finite source at its parent. Its selected value stays
unresolved until recursive materialisation or realisation reaches it. Nested
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

## Eager public construction

Gridspace delays selection, not ordinary construction. Each public builder
applies one rule:

```text
no direct Grid or Gridspace input  -> construct the eager value now
at least one explicit source       -> preserve sources, singleton-wrap siblings,
                                      and return a Gridspace
```

The current behaviour is:

| Entry point | Scalar or complete input | Explicit varying input |
|---|---|---|
| `Material` | `Materials.Material` | `Gridspace{Materials.Material}` |
| cable-part declaration | callable `PartBuilder` | `Gridspace{PartBuilder}` |
| `CableBuilder` | `DataModel.CableDesign` | `Gridspace{DataModel.CableDesign}` |
| `at`, `trifoil`, `hflat`, `vflat` | `PositionDefinition` | `Gridspace{PositionDefinition}` |
| `Earth` | `EarthModel` | `Gridspace{EarthModel}` |
| `SystemBuilder` | `LineParametersProblem` | `Gridspace{LineParametersProblem}` |
| `CableConstantsProblem` | `CableConstantsProblem` | `Gridspace{CableConstantsProblem}` |
| `@gridspace` keyword constructor | strict struct | `Gridspace{Target}` |

A cable part remains callable because its absolute radius depends on the
preceding radial layers. The callable constructs one radial part. The
other scalar-complete builders construct their domain values immediately.

## Callable builders

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
for local experiments. Exported builders use concrete callable types.
Do not introduce a passive record merely to store arguments that another
function immediately unpacks.

`@gridspace` applies the same rule to a keyword-constructed struct. The macro
retains the strict positional constructor, returns the struct immediately for
scalar keyword input, and creates a Gridspace only when a field is an explicit
finite source. `@gridspace` composes with `@relax` in either order.

## Materialisation and realisation

Ordinary iteration selects an internal unresolved point and recursively
materialises its arguments. Deterministic values pass through unchanged.
nested points invoke their own callable before the parent callable is invoked.
After loading Measurements, an `UncertainValue` materialises as one
`Measurement`.

Stochastic realisation follows the same recursion with a caller-owned random
number generator. Only `UncertainValue` leaves are redrawn. Deterministic
selections remain fixed. A zero-sigma descriptor resolves deterministically.

Materialisation and realisation are internal. Developers extend the public grammar through
concrete builders and supported uncertainty extensions, not
by exposing unresolved points as application data.

## Higher-order computation

Combinatorial traversal has one sequence:

```text
select point -> `materialize` complete core problem
             -> Engine.compute
             -> append core result
```

`ParametricResult` stores the higher-order formulation and the ordered vector
of core results. The result does not retain traversal state.

Direct linear propagation uses the same traversal. The Measurements extension
changes only how an uncertain descriptor materialises. `LinearErrorResult`
likewise stores only its formulation and ordered core results.

Monte Carlo selects each outer point once, derives a deterministic point seed,
and repeatedly realises that same point:

```text
for each selected outer point
    for each trial
        redraw uncertain leaves within that point
        build a fresh complete core problem
        Engine.compute
    aggregate that point's draws
end
```

Multiple nominal/error points therefore produce multiple aggregates, not one
mixture. `MonteCarloResult` directly owns sample-mean core results, statistics,
optional retained samples, optional histograms, the root seed, point seeds,
and trial counts.

All completed result spaces are one-dimensional finite Julia collections.
Iteration and indexing return one stored core result per original Gridspace
point, in traversal order. Monte Carlo iteration returns the representative
core result reconstructed from each point's sample means; individual trials
remain available only through `samples`. Standard `first`, `last`, `only`,
`collect`, `map`, and `zip` operations apply. `only` asserts singleton
cardinality and performs no statistical selection or projection.

## Projecting completed result spaces

[`project`](@ref) converts a completed result space into the next finite space
of complete problems. The projection definition owns the mathematical choice;
ParametricBuilder owns the fixed sequence:

```text
entitle -> select -> derive -> materialize -> finish
```

`entitle` rejects an unsupported definition/result-space pair before partial
work. `select` aligns the result products required by the projection. `derive`
chooses or calculates a finite set of representative states. `materialize`
constructs one complete downstream problem per representative. The derived
`finish` method checks that the output is nonempty and type-consistent, then
returns `Gridspace{P}` in representative order.

Only `project` is exported from the package root. An external implementation
imports the definition and required stages from their owner:

```julia
import LineCableModels.ParametricBuilder:
    AbstractProjectionDefinition,
    entitle,
    select,
    derive,
    materialize

struct EnvelopeProjection <: AbstractProjectionDefinition end

entitle(::EnvelopeProjection, source::MonteCarloResult) = source
select(::EnvelopeProjection, source::MonteCarloResult) = (
    results = source,
    statistics = statistics(source),
)
# Define derive and materialize for the extension's representative and problem types.
```

A summary projection may calculate a synthetic representative whose components
did not occur together in any trial. A retained-trial projection instead
selects complete joint samples and preserves their correlation. The definition
must state which operation it performs. Projection never overloads `only` and
does not add a result wrapper; its completed product is the next `Gridspace`.

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

## Optional package extensions

The core package declares uncertainty without loading Measurements or
Distributions.

- Loading Measurements adds direct materialisation of `UncertainValue` while
  retaining exact structural reuse.
- Loading Distributions adds standardised univariate sampling families. The
  selected distribution must have finite mean and positive finite standard
  deviation. Samples are transformed to the descriptor's nominal value and
  standard uncertainty.

Neither extension knows about finite-source identity, result presentation, or
Engine internals.

## Performance and conformance

The implementation relies on tuple-specialised recursion and Julia's public
product and zip iterators. The implementation guarantees:

- `length` is analytic and never enumerates points.
- Product and zip traversal are linear in yielded work.
- Materialisation and realisation use no dictionary or identity lookup.
- Immutable scalar targets infer through selection, materialisation, and
  realisation.
- After warmup, deterministic iteration and a bare 10,000-realisation scalar
  loop add zero heap allocations.
- Full cable and line construction may allocate only what their existing
  vectors, domain constructors, Engine computation, and requested result
  storage intrinsically require.

The conformance suite in `test/unit/parametricbuilder/conformance.jl` checks
these properties, exact structural reuse, explicit variation, eager public
construction, and the absence of a random-access Gridspace API.

## Implementation map

The implementation is split across:

- `src/parametricbuilder/grid.jl`: finite values and uncertainty descriptors.
- `src/parametricbuilder/gridspace.jl`: composition, point selection, recursive
  materialisation, and realisation.
- `src/parametricbuilder/macros.jl`: strict scalar construction and explicit
  Gridspace lifting for `@gridspace`.
- material, cable, position, and system builder files: eager construction
  rules and concrete callable algorithms.
- `src/parametricbuilder/compute.jl`: combinatorial traversal.
- `src/uq/linearerror.jl` and `src/uq/montecarlo/compute.jl`: direct and
  repeated stochastic traversal.
- Measurements and Distributions extensions: dependency-specific uncertainty
  behaviour only.
