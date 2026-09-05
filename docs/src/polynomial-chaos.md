# Polynomial chaos

`PolynomialChaos` performs non-intrusive uncertainty propagation: it selects
physical input values, calls the unchanged deterministic cable solver at a
tensor quadrature rule, fits native output expansions, and checks every fitted
surrogate against independent core solves. The core package defines the public
formulation and result contracts without loading PolyChaos.jl.

## Activation and construction

Install and load PolyChaos alongside LineCableModels:

```julia
using LineCableModels
using PolyChaos
```

Loading PolyChaos activates `LineCableModelsPolyChaosExt`. The public call is:

```julia
result = compute(
    ParametricProblem(problem_space),
    PolynomialChaos(
        CableConstantsFormulation();
        degree = 3,
        quadrature_order = nothing,
        distribution = :normal,
        validation_points = 16,
        validation_rtol = 1e-3,
        max_evaluations = 10_000,
        validation_seed = 0,
        options = (retain_details = false,),
    ),
)
```

`quadrature_order=nothing` resolves to `degree + 1`. Without `using
PolyChaos`, the types remain available for documentation and configuration,
but the `compute(::ParametricProblem, ::PolynomialChaos)` method is absent.
There is no Monte Carlo fallback.

## Input probability model

Every positive-uncertainty descriptor is one independent stochastic
coordinate. Exact structural reuse of one builder argument remains one
coordinate even when that argument is used several times. Distinct descriptor
slots remain independent even when their values are equal. Correlation and
covariance matrices are not accepted in this implementation.

| `distribution` | latent measure | physical interpretation |
|:--|:--|:--|
| `:normal` | ``\xi \sim \mathcal{N}(0,1)`` | declared nominal and uncertainty are mean and standard deviation |
| `:uniform` | ``\xi \sim \mathcal{U}(-1,1)`` | an affine uniform variable with the declared mean and standard deviation |

A zero-uncertainty descriptor keeps its positional slot and is supplied at its
nominal value, but it does not add a stochastic coordinate. A point with no
active coordinates is rejected before a core solve; use `Combinatorial` for a
deterministic space. Physical values are never clamped or repaired. If the
declared distribution reaches invalid radii, resistivities, or other inputs,
the owning constructor or validator throws normally.

## Cost and validation

For active dimension ``d``, total degree ``p``, quadrature order ``q``, and
``V`` validation points, one outer point has

```math
M = \binom{d+p}{p}, \qquad Q = q^d, \qquad N_{solve} = Q + V.
```

`M` is the number of total-degree basis terms and `Q` is the complete tensor
collocation count. Before the first core solve, the extension computes these
costs for every outer Gridspace point and rejects the whole request when their
sum exceeds `max_evaluations`. It never reduces the degree/order, switches
method, or continues after a warning.

Validation is mandatory. A local RNG derived from `validation_seed` and the
outer-point index generates independent latent points from the selected
measure. The extension compares predicted and actual native quantities with
relative RMS and relative maximum errors. Both errors must be at most
`validation_rtol` for every quantity. A failure reports the point, stochastic
dimension, degree/order, basis and solve counts, worst quantity and native
index/frequency, measured errors, and requested tolerance. The global RNG is
not used.

This full tensor method is intended for low-dimensional, smooth responses. It
is not universally faster than Monte Carlo: ``q^d`` grows exponentially, and
nonsmooth or poorly resolved responses are rejected by the budget or
validation guard.

## Results and scientific observations

`PolynomialChaosResult` is aligned with outer Gridspace order. Ordinary
iteration and indexing return mean `CableConstants` or `LineParameters` core
results. Product accessors return aligned fitted data:

```julia
using Statistics

mean_core = only(result)
moments = only(statistics(result))
expansion = only(expansions(result))
diagnostic = only(validation(result))
retained = details(result)

mean_R = observe(result, statistics, R, Statistics.mean, 1, 1)
std_R = observe(result, statistics, R, Statistics.std, 1, 1)
```

Each expansion stores its PolyChaos basis and one fixed R/L/C/G coefficient
product. The polynomial-term axis is last. Each statistics entry stores only
native `(mean, std)` values; it is deliberately not a `SampleSummary`.
`observables` applies display units only when publishing explicit mean or
standard-deviation requests.

The solver fits R/L/C/G directly. For line results it reconstructs
``Z=R+j\omega L`` and ``Y=G+j\omega C`` while preserving domain, frequency
axis, matrix shape, ordering, and `:pul`/`:total` basis. Fitting Z/Y magnitude
or angle would mix coordinate choices, introduce phase wrapping, and obscure
the physical quantities used by the existing observation API.

`samples(result)` and `histograms(result)` are intentionally undefined. Cheap
values generated later with `PolyChaos.evaluatePCE` are surrogate evaluations,
not retained deterministic core-solver trials and not empirical histograms.
The non-CI benchmark in `benchmarks/polynomial_chaos.jl` demonstrates such
offline evaluation while separately reporting actual core solve counts. Its
Monte Carlo reference defaults to 5,000 trials; set
`LINECABLEMODELS_PCE_BENCHMARK_TRIALS` to change that workload.

When `retain_details=true`, `details(result)` has the fixed schema
`(points=[(collocation=[...], validation=[...]), ...],)`, containing the core
computation details for every actual solve. Otherwise it is the empty named
tuple.
