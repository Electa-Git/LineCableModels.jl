# Global sensitivity

The optional GlobalSensitivity extension calculates Sobol first-order,
total-order, and optional second-order indices for native LineCableModels
observable requests. It repeatedly realises the uncertainty descriptors in
each selected Gridspace point, evaluates the ordinary core formulation, and
attributes output variance to the active input coordinates.

Install and load the three optional dependencies together:

```julia-repl
pkg> add GlobalSensitivity QuasiMonteCarlo Distributions
```

```julia
using LineCableModels
using Distributions
using GlobalSensitivity
using QuasiMonteCarlo
```

The packages activate an extension; loading `LineCableModels` alone does not
load them.

## Workflow

Start with a `Gridspace{LineParametersProblem}` or other supported parametric
problem space containing `UncertainValue` descriptors. Select outputs with the
same native requests used by `observe` and `observables`:

```julia
problem = ParametricProblem(problem_space)

analysis = Sensitivity(
    Formulation(),
    GlobalSensitivity.Sobol(),
    (
        @observe(R[1, 1, 1]),
        @observe(L[:, :, 1]),
    );
    samples=4096,
    distribution=:normal,
    input_labels=("core radius", "insulation thickness"),
    max_evaluations=16_384,
    max_output_values=81_920,
)

sensitivity = compute(problem, analysis)
product = only(sensitivity)

product.inputs
product.active_indices
product.evaluations
first_order(sensitivity)
total_order(sensitivity)
```

`Sensitivity` preserves the supplied `GlobalSensitivity.Sobol` method and
QuasiMonteCarlo sampler. The default sampler is
`QuasiMonteCarlo.SobolSample()` and the default first-order estimator is
`:Jansen1999`. The other admitted estimators are `:Homma1996`, `:Sobol2007`,
and `:Janon2014`.

Use `GlobalSensitivity.Sobol(order=[0, 1, 2])` to request second-order
products. The supported order sets are `[0, 1]` and `[0, 1, 2]`. With
`nboot=1`, confidence products are `nothing`; larger bootstrap counts return
arrays aligned with the corresponding indices.

## Physical input mapping

Each nonzero-standard-uncertainty descriptor is one active coordinate.
`:normal` maps its unit-cube coordinate through a normal distribution with the
declared nominal value and standard uncertainty. `:uniform` uses the bounded
interval `nominal ± sqrt(3) * sigma`, which has that same standard
uncertainty. Zero-sigma descriptors retain their nominal value and are omitted
from the sensitivity dimension.

`input_labels` must cover every uncertainty descriptor before zero-sigma
coordinates are removed. A result stores the labels of active coordinates in
`inputs` and their original one-based descriptor positions in
`active_indices`. Without explicit labels, active coordinates are named `x1`,
`x2`, and so on.

Distinct descriptor positions are independent. Exact structural reuse within
one builder argument remains one coordinate. General correlated inputs require
a dependency-owned joint sensitivity method; the Sobol integration does not
infer covariance from equal values or reused `Grid` objects.

## Products and shapes

A [`SensitivityResult`](@ref) is a one-dimensional collection with one product
per outer Gridspace point. Within each product, `first_order`, `total_order`,
and `second_order` are tuples in request declaration order.

For a request whose observed value has shape `output_shape`, first- and
total-order products have shape `(output_shape..., d)`. Second-order products
have shape `(output_shape..., d, d)` and store each input pair once in the
upper triangle; the lower triangle is zero. A scalar request therefore returns
a length-`d` vector for first and total order, and a `d × d` matrix for second
order.

The point-aligned accessors are:

```julia
first_order(sensitivity)
total_order(sensitivity)
second_order(sensitivity)
confidence(sensitivity)
```

Setting `options=(retain_details=true,)` retains one core computation-details
record per evaluation under `details(sensitivity).evaluations`. The records
remain aligned first by outer point and then by core evaluation.

## Cost and failure boundaries

Let `N` be `samples`, `B` be `method.nboot`, and `d` be the number of active
coordinates at one outer point. First- and total-order analysis performs

```text
N * B * (d + 2)
```

core evaluations. Enabling second order performs

```text
N * B * (2d + 2)
```

evaluations. The extension calculates the cost of every outer point before the
first solve and rejects totals above `max_evaluations`. It probes one result to
establish request shapes, then rejects flattened storage above
`max_output_values`.

Every observed value must retain its type and shape and contain finite real
numbers. Empty, complex, nonfinite, constant, or shape-changing observations
fail explicitly. Wrapped angle observations are not admitted for variance
decomposition. Physical construction errors, including `DomainError`, are
propagated; Sobol analysis does not reject and redraw infeasible points.

QuasiMonteCarlo 0.3.6 emits its documented warning when the deterministic
`SobolSample()` design-matrix overload uses `NoRand()`. The integration keeps
that public deterministic API and does not hide the warning or substitute a
private randomisation path.

## Analytical benchmark

The non-CI Ishigami benchmark checks first-, total-, and second-order products
against their closed-form reference and reports the exact evaluation count,
elapsed time, and maximum errors:

```bash
julia --project=test benchmark/sobol_ishigami.jl
```

The benchmark is executable maintenance evidence. It is not discovered by the
package test runner and is not a substitute for the deterministic analytical
and real cable-domain test suites.
