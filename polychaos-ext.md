# LineCableModels Non-Intrusive Polynomial Chaos Integration

You are working in **implementation mode** on the live `LineCableModels.jl` repository.

Implement a non-intrusive polynomial-chaos UQ formulation backed by `PolyChaos.jl`. The formulation must evaluate the existing deterministic LineCableModels core solver only at the required stochastic collocation nodes, fit polynomial-chaos expansions of the native scientific outputs, validate those expansions against independent core solves, and return a proper higher-order result aligned with the source `Gridspace`.

This is a goal-pursuit task. Do not return another architecture essay or planning document. Inspect the live repository, reconcile this contract against current paths and already-completed refactors, execute the finite commit sequence, run the repository-native validation gates, and stop only on an explicit stop condition.

---

## 1. Authority and preflight

Before modifying source:

1. Read the live repository instructions, maintained architecture/conventions documentation, current `README*`, UQ documentation, `Project.toml`, CI, test runner, and current `git status`.
2. Read the current implementations of:
   - `Grammar.compute`, result-space abstract types, `observe`, `observables`, `ComputationOptions`, and `ComputationDetails`;
   - `ParametricBuilder.Grid`, `UncertainValue`, `GridPoint`, `Gridspace`, `points`, `uncertainties`, `materialize`, and `realize`;
   - `Combinatorial`, `MonteCarlo`, `LinearError`, and their result containers;
   - `CableConstants`, `LineParameters`, and the native `R/L/C/G/Z/Y` observation methods;
   - optional dependency extensions and extension tests.
3. Inspect recent relevant commits. Treat the live branch as authoritative when prior snapshots differ.
4. Run the repository's current green baseline before editing.
5. Record unrelated modified/untracked files and preserve them exactly.
6. Inspect the **resolved PolyChaos release and source** used by the project. The intended dependency line is PolyChaos `2.1`; do not copy undocumented internals or signatures from stale examples.

Resolve conflicts in this order:

1. this prompt and explicit current user requirements;
2. Julia Core/Base and documented dependency interfaces;
3. the repository's Native-First Semantic Economy policy;
4. ownership-centered recursive module layout;
5. dispatch-driven Template Method policy;
6. current maintained repository documentation and supported tests;
7. implementation details and historical artifacts.

Do not overwrite or reverse unrelated in-progress result-observation, ReportBuilder, PlotBuilder, or result-space reconciliations. Integrate with their **live accepted contracts**.

---

## 2. Preservation contract

Preserve all existing behavior not explicitly changed here:

- `Grid`, `Gridspace`, `GridPoint`, product and zip semantics;
- outer-Gridspace traversal order;
- deterministic `materialize(point)` behavior, including Measurements-based propagation;
- existing `Combinatorial`, `LinearError`, and `MonteCarlo` public calls and result values;
- Monte Carlo distributions, seeds, point-seed derivation, trial counts, retained samples, histograms, details, and DKW sizing;
- core CableConstants and LineParameters numerical algorithms and baselines;
- `observe`, `observables`, `Units`, `:pul`/`:total`, plotting, and report behavior;
- current result-space collection behavior and core-result `eltype` semantics;
- current public API not explicitly extended below.

The existing materialisation/realisation contract is authoritative and must remain intact. For the same root seed, Gridspace, distribution, and Julia/dependency versions, Monte Carlo must preserve its current random-draw order, results, seeds, trial counts, retained products, and failure behavior. This task extends `realize`; it does not replace, rename, or bypass it.

---

## 3. Locked architectural vocabulary

The repository already owns two distinct internal actions. Preserve their meanings exactly:

```text
materialize
    recursively construct the deterministic selected point;
    deterministic leaves pass through;
    nested GridPoints invoke their callable;
    UncertainValue materialises through Measurements when that extension is loaded

realize
    recursively construct one concrete stochastic/collocation instance;
    deterministic leaves remain fixed;
    only UncertainValue leaves are resolved;
    the caller owns the random generator or supplies the concrete physical values
```

The code vocabulary is:

```text
Grid
GridPoint
Gridspace
points
uncertainties
materialize
realize
compute
PolynomialChaos
PolynomialChaosResult
statistics
expansions
validation
details
observe
observables
```

The public higher-order action remains:

```julia
compute(problem, formulation)
```

The stable parameter-space sequence for polynomial chaos is:

```text
points(space)
    -> GridPoint
    -> uncertainties(point)
    -> PolyChaos-owned latent nodes and weights
    -> formulation-owned mapping to concrete physical values
    -> realize(point, values)
    -> complete core problem
    -> compute(core_problem, inner)
    -> formulation-owned spectral projection and validation
```

`materialize(point)` remains the deterministic/Measurements path. `realize(rng, point, distribution)` remains the Monte Carlo path. `realize(point, values)` is the new deterministic-node overload used by polynomial chaos and later supported methods such as quasi-Monte Carlo or sparse-grid collocation.

Materialisation and realisation remain internal ParametricBuilder protocols. Do not package-root export them and do not expose `GridPoint` or unresolved points as application data.

Do not introduce any of the following:

```text
Realization
AbstractRealization
RealizationPolicy
RandomRealization
CollocationRealization
RealizationMapper
MaterializationPolicy
```

The existing verb `realize` is the semantic authority. Extend it with methods; do not create a hierarchy around it.

Do not introduce a LineCableModels-owned generic named:

```text
sample
generate
nodes
collocate
parameterize
propagate
run_uq
```

Random generation remains owned by `rand`; PolyChaos basis, node, weight, and PCE evaluation remain owned by PolyChaos. `realize` owns recursive replacement of uncertainty leaves with concrete values. `materialize` owns deterministic Measurements-based construction. `compute` owns execution.

---

## 4. ParametricBuilder realisation extension

This task does **not** redesign Gridspace. It admits one missing structural query and one method on the existing `realize` generic so deterministic collocation methods can use the same recursion without impersonating Monte Carlo.

### 4.1 Internal owner-visible protocol

Within `LineCableModels.ParametricBuilder`, retain the existing non-root-exported protocol:

```julia
GridPoint
points
materialize
realize
```

Add:

```julia
uncertainties
```

These names are internal owner-visible extension seams. `UQ` and optional package extensions import or qualify them explicitly from `LineCableModels.ParametricBuilder`. Do not add them to the package-root export list and do not advertise unresolved points as ordinary user data.

Preserve all existing methods and semantics, including:

```julia
materialize(value)
materialize(value::UncertainValue)
materialize(point::GridPoint)

realize(rng::Random.AbstractRNG, value, distribution)
realize(rng::Random.AbstractRNG, value::UncertainValue, distribution)
realize(rng::Random.AbstractRNG, point::GridPoint, distribution)
```

The Measurements extension remains the owner of deterministic `UncertainValue -> Measurement` materialisation. The Distributions extension/current Monte Carlo path remains the owner of random `UncertainValue` draws.

### 4.2 `uncertainties`

Add one admitted generic:

```julia
uncertainties(point::GridPoint)
```

Its exact meaning is:

> Return a concrete tuple containing every `UncertainValue` leaf of one selected `GridPoint`, in the same deterministic depth-first, left-to-right order used by recursive realisation.

Required behavior:

```julia
uncertainties(::Any) = ()
uncertainties(value::UncertainValue) = (value,)
uncertainties(point::GridPoint) = ...
```

The `GridPoint` method recursively traverses nested `GridPoint`s and concatenates tuples. It does **not** descend into arbitrary tuple- or array-valued deterministic leaves, because existing materialisation/realisation does not do so. Do not return a `Vector`, dictionary, identity key, binding record, coordinate record, or configuration object.

Structural laws:

- traversal order is deterministic;
- the order matches `realize(point, values)` exactly;
- nested GridPoints use the same recursion;
- one uncertainty occupying one builder argument occupies one coordinate, even when the builder reuses that argument internally;
- two uncertainty slots occupy two coordinates even if their descriptor values are equal or they originated from the same Grid object;
- no object-identity deduplication occurs;
- zero-sigma descriptors remain in the tuple;
- changing unrelated deterministic values does not reorder uncertainty coordinates.

This generic is admitted because the coordinate order is real ParametricBuilder semantics required by more than one uncertainty method. Neither Base nor PolyChaos can infer it from LineCableModels builders.

### 4.3 `realize(point, values)`

Add a method to the existing generic:

```julia
realize(point::GridPoint, values::Tuple)
```

Its exact meaning is:

> Recursively construct the complete target owned by `point`, replacing every `UncertainValue` leaf with the corresponding concrete physical value from `values`.

Requirements:

- verify `length(values) == length(uncertainties(point))` before invoking any builder;
- consume values in the exact `uncertainties(point)` order;
- recurse through nested GridPoints;
- deterministic leaves pass through unchanged;
- zero-sigma descriptors are supplied as their nominal values by the owning formulation;
- invoke each nested callable before its parent callable, exactly like existing materialisation and random realisation;
- return an instance of the same Gridspace target family; concrete numeric parameterization may differ according to the supplied values, exactly as current materialisation and random realisation already permit;
- let owned constructors/validation reject invalid physical values;
- do not clamp, coerce, repair, or silently substitute values;
- keep the recursive implementation owner-local; use local recursive functions unless inference/allocation evidence warrants one narrowly named internal kernel;
- do not add a cursor, context, mapper, realization object, or result wrapper.

For an immutable scalar target after compilation, representative calls to `uncertainties(point)` and `realize(point, values)` must infer successfully and add zero heap allocations beyond the target constructor itself.

### 4.4 Preserve Monte Carlo unchanged

Do not migrate Monte Carlo merely to make it use the new tuple overload. Its current method:

```julia
realize(rng, point, distribution)
```

is already the correct stochastic action and may retain its direct recursive implementation for hot-path performance and exact draw ordering.

Required preservation:

- same draw order for the same seed and point structure;
- same zero-sigma behavior;
- same result types and numerical values;
- same point seeds, root seed, trial counts, samples, histograms, statistics, and details;
- same error behavior.

Add tests proving that `realize(point, physical_values)` and `realize(rng, point, distribution)` follow the same structural leaf order, without rewriting the random path through an allocating temporary tuple.

---

## 5. Optional PolyChaos dependency boundary

Add PolyChaos as an optional dependency:

```toml
[weakdeps]
PolyChaos = "8d666b04-775d-5f6e-b778-5ac7c70f65a3"

[extensions]
LineCableModelsPolyChaosExt = "PolyChaos"

[compat]
PolyChaos = "2.1"
```

Add PolyChaos to test extras/targets using the repository's current extension-test convention.

Create the optional extension under the ownership-centered layout:

```text
ext/LineCableModelsPolyChaosExt/
├── LineCableModelsPolyChaosExt.jl
├── basis.jl
├── projection.jl
└── compute.jl
```

Merge files when a responsibility remains too small to deserve a separate file. Do not create `helpers.jl`, `utils.jl`, `common.jl`, `adapter.jl`, or another catch-all file.

Core source must not import PolyChaos. Core formulation/result types remain parametric over extension-owned basis and coefficient payload types.

When PolyChaos is not loaded:

- LineCableModels must load and precompile normally;
- `PolynomialChaos` and `PolynomialChaosResult` remain defined and documentable;
- no fake fallback computation is provided;
- unsupported `compute(problem, PolynomialChaos(...))` fails through ordinary missing extension dispatch rather than silently falling back to Monte Carlo.

---

## 6. Public formulation contract

Define in `UQ`, without importing PolyChaos:

```julia
struct PolynomialChaos{F, O <: ComputationOptions} <: AbstractFormulation
    inner::F
    degree::Int
    quadrature_order::Int
    distribution::Symbol
    validation_points::Int
    validation_rtol::Float64
    max_evaluations::Int
    validation_seed::UInt64
    options::O
end
```

Provide the public constructor:

```julia
PolynomialChaos(
    inner;
    degree = 3,
    quadrature_order = nothing,
    distribution = :normal,
    validation_points = 16,
    validation_rtol = 1e-3,
    max_evaluations = 10_000,
    validation_seed = 0,
    options = (;),
)
```

Normalize `quadrature_order = nothing` to `degree + 1` before constructing the object.

Constructor laws:

- `degree >= 1`;
- `quadrature_order >= degree + 1`;
- `distribution in (:normal, :uniform)`;
- `validation_points >= 1`;
- `validation_rtol > 0` and finite;
- `max_evaluations >= 1`;
- `validation_seed` converts exactly to `UInt64`;
- `options` is normalized through:
  ```julia
  computation_options(Val(PolynomialChaos), options)
  ```
- the initial computation-owned option set is exactly:
  ```julia
  (retain_details = false,)
  ```
- unknown options fail immediately.

Do not add a `sampler`, `basis`, `method`, `mode`, covariance matrix, adaptive-degree flag, fallback formulation, or generic options dictionary in this first implementation.

Export `PolynomialChaos` from `LineCableModels.UQ` and the package root.

---

## 7. Result contract

Define in `UQ`, without importing PolyChaos, one result-space type following the live repository's current result hierarchy and Base collection interface:

```julia
struct PolynomialChaosResult{
    T,
    F,
    E <: AbstractVector,
    S <: AbstractVector,
    V <: AbstractVector,
    D <: ComputationDetails,
} <: AbstractUncertaintyResult{T}
    formulation::F
    values::Vector{T}
    expansion_values::E
    stats::S
    validation_values::V
    details::D
end
```

Adapt the supertype only if the live branch has already established `AbstractResultSpace`; do not create a parallel hierarchy.

Constructor invariants:

- `values`, expansions, statistics, and validation records contain one entry per outer `GridPoint`;
- supplemental details are empty or follow one fixed owner-defined schema;
- nested result-space envelopes remain rejected according to the live core-result guardrail;
- result element type is exposed through standard `eltype` and ordinary iteration;
- do not overload `Base.only`; Base must work through the collection interface.

Add and export these UQ product accessors:

```julia
expansions(result::PolynomialChaosResult)
validation(result::PolynomialChaosResult)
```

Reuse the existing:

```julia
statistics(result::PolynomialChaosResult)
details(result::PolynomialChaosResult)
```

Each `expansions(result)[point]` is one fixed-key concrete named tuple:

```julia
(
    basis = polychaos_basis,
    coefficients = coefficient_product,
)
```

Each coefficient product is a concrete named tuple:

```julia
# CableConstants
(R = ..., L = ..., C = ...)

# LineParameters
(R = ..., L = ..., C = ..., G = ...)
```

The polynomial-term axis is the final array dimension. Scalar cable constants use coefficient vectors.

Each `statistics(result)[point]` is a concrete named tuple with the same quantity keys. Each value contains exactly native mean and standard-deviation products, for example:

```julia
R = (mean = ..., std = ...)
```

Do not reuse `SampleSummary`, because polynomial chaos does not directly own empirical minima, quantiles, maxima, or trial counts.

Each validation record is a fixed-key concrete named tuple containing at least:

```text
stochastic dimension
polynomial degree
basis term count
collocation evaluation count
validation evaluation count
per-observable relative RMS error
per-observable relative maximum error
worst observable/index/frequency location
```

Do not use `Dict{Symbol,Any}`.

Do not overload `samples` or `histograms` for `PolynomialChaosResult`. Surrogate-generated values are not retained core-solver trials and are not empirical histograms.

---

## 8. Supported initial domain

The initial PolyChaos extension supports only outer Gridspaces whose core calculations return:

```text
CableConstants
LineParameters
```

The extension may use one owner-local, non-exported algorithmic generic:

```julia
fit_expansion(..., results::AbstractVector{<:CableConstants})
fit_expansion(..., results::AbstractVector{<:LineParameters})
```

This generic is explicitly admitted under the algorithm and dispatch-contract warrants:

- it owns non-intrusive spectral projection of one core-result family;
- two current core-result types require materially different extraction, invariants, reconstruction, and coefficient shapes;
- it remains inside `LineCableModelsPolyChaosExt`;
- it is not package-root exported;
- no sibling or unrelated extension imports an underscore-prefixed substitute.

Do not add separate module-level helpers for:

```text
extract_outputs
flatten_outputs
rebuild_result
mean_result
fit_quantity
```

Keep repeated linear algebra direct or use local functions inside each `fit_expansion` method.

Unsupported core result families fail through ordinary dispatch at the PolyChaos extension boundary. Do not add a broad fallback returning `nothing`, `Any`, or an empty expansion.

---

## 9. PolyChaos dependency contract

Use only public PolyChaos `2.1` interfaces verified from the resolved source.

The implementation must account for these actual API facts:

1. `GaussOrthoPoly(degree; Nrec = ...)` and `Uniform_11OrthoPoly(degree; Nrec = ...)` attach a Gaussian quadrature rule with `Nrec - 1` points.
2. To obtain `quadrature_order` points, construct each univariate basis with:
   ```julia
   Nrec = quadrature_order + 1
   ```
3. `MultiOrthoPoly(univariate_bases, degree)` creates the total-degree multivariate basis.
4. `PolyChaos.nw(mop)` returns one node vector and one weight vector **per stochastic coordinate**. It does not return the full tensor-product node matrix.
5. Build the complete tensor rule explicitly and deterministically with `Iterators.product`, using Julia product order.
6. `PolyChaos.evaluate(node_matrix, mop)` expects a matrix shaped:
   ```text
   number_of_nodes × stochastic_dimension
   ```
   and returns:
   ```text
   number_of_basis_terms × number_of_nodes
   ```
7. `PolyChaos.convert2affinePCE` supports conversion from mean/standard deviation for the standard Gaussian and `Uniform_11` bases.
8. `PolyChaos.evaluatePCE` supports multivariate evaluation of one coefficient vector at a matrix of latent points.

Do not access PolyChaos private fields or private functions when a public method exists. Do not add methods to undocumented PolyChaos internals.

---

## 10. Stochastic coordinate model

For every selected outer `GridPoint`:

```julia
descriptors = uncertainties(point)
```

Require every descriptor's nominal value and standard uncertainty to be finite real scalars.

A descriptor with zero standard uncertainty:

- remains in the positional `uncertainties(point)` contract;
- does not create a stochastic coordinate;
- is supplied to `realize(point, values)` at its nominal value.

The active stochastic dimension is the number of descriptors with positive standard uncertainty.

Reject a point with zero active stochastic dimension before any core solve and instruct the caller to use `Combinatorial` for a deterministic space.

Initial probability support is exactly:

### Normal

```text
latent coordinate ξ ~ Normal(0, 1)
physical value x has declared mean μ and standard deviation σ
basis: GaussOrthoPoly
```

### Uniform

```text
latent coordinate ξ ~ Uniform(-1, 1)
physical value x has declared mean μ and standard deviation σ
basis: Uniform_11OrthoPoly
```

Use `PolyChaos.convert2affinePCE` to obtain the affine mean/standard-deviation mapping for each active descriptor. Evaluate that affine mapping at each latent coordinate using the public PolyChaos basis/evaluation contract or an exactly equivalent allocation-free affine evaluation verified against `evaluatePCE`.

Do not clamp, truncate, reflect, or repair invalid physical values produced by a declared distribution. If a Gaussian tail creates invalid resistivity, radius, permittivity, or another domain input, the declared stochastic model is invalid and the core constructor/validator must reject it.

Initial multivariate semantics are:

```text
independent active coordinates
+ exact structural reuse when one builder argument is used in several places
```

Do not accept a covariance/correlation matrix in this campaign. PolyChaos `MultiOrthoPoly` represents a product measure. General correlated inputs require a separately designed independent-latent transform and are outside this slice.

---

## 11. Cost preflight

For every point, let:

```text
d = active stochastic dimension
p = polynomial degree
q = quadrature order
M = binomial(d + p, p)          # total-degree basis terms
Q = q^d                         # tensor collocation evaluations
V = validation_points
```

Before executing any core solve:

1. collect or traverse every outer `GridPoint` without materializing core problems;
2. calculate `d`, `M`, `Q`, and `V` for every point;
3. calculate the total requested core evaluations:
   ```text
   sum(Q + V for every outer point)
   ```
4. throw an `ArgumentError` when the total exceeds `max_evaluations`;
5. include outer-point count, per-point dimensions, basis terms, collocation evaluations, validation evaluations, requested total, and budget in the error.

The budget is a hard stop. Do not silently:

- lower the degree;
- lower quadrature order;
- reduce validation points;
- switch to regression;
- switch to sparse grids;
- switch to Monte Carlo;
- continue after warning.

No core compute call may occur before this preflight succeeds.

---

## 12. Collocation evaluation

For each outer point:

1. build one univariate basis for each active descriptor;
2. build one `MultiOrthoPoly` total-degree basis;
3. obtain per-coordinate nodes/weights through `PolyChaos.nw`;
4. form one complete tensor-product latent-node matrix and one product-weight vector;
5. verify:
   - node matrix shape is `Q × d`;
   - weight length is `Q`;
   - all nodes and weights are finite;
   - all weights are positive;
   - product weights sum to one within a strict numerical tolerance appropriate to `Float64`;
6. construct each complete physical-value tuple in the full `uncertainties(point)` order;
7. call:
   ```julia
   core_problem = realize(point, physical_values)
   core_result = compute(
       core_problem,
       formulation.inner;
       options = problem.options,
   )
   ```
8. require all collocation results for one point to share one concrete result type and compatible result metadata.

Do not create a second Gridspace by taking a Cartesian product of the PolyChaos node columns. The tensor node rows are already complete stochastic coordinates. Accidentally applying a second product is an execution-budget catastrophe and must be guarded by tests.

Retained computation details, when requested, are captured through the live `ComputationDetails` owner contract. Do not import a private underscore-prefixed computation-owner helper across modules. Use the live owner-visible contract, or the direct Base type-token expression already sanctioned by the repository.

---

## 13. Spectral projection

For one point, let:

```text
Φ = PolyChaos.evaluate(nodes, basis)     # M × Q
w = product quadrature weights          # Q
Y = flattened native core outputs       # Q × N
```

Compute basis norms using the same quadrature rule:

```text
hα = Σq wq Φ[α,q]^2
```

Reject nonfinite or nonpositive norms.

Fit all scalar outputs of one quantity family in one matrix operation:

```text
C = H⁻¹ Φ W Y
```

with shapes:

```text
C: M × N
H: diagonal M × M, implemented without constructing a dense diagonal matrix
W: diagonal Q × Q, implemented without constructing a dense diagonal matrix
```

Do not fit one scalar output at a time. Matrix entries and frequencies are flattened into output columns, projected together, and reshaped afterward.

Find the constant basis term from the unique all-zero row of the PolyChaos multi-index matrix. Do not assume its position without asserting the invariant.

For an orthogonal probability basis:

```text
mean = constant coefficient
variance = Σ(nonconstant α) Cα² hα
std = sqrt(max(variance, 0))
```

The zero clamp is permitted only for negative variance on the scale of floating-point roundoff. It is not permission to rewrite physical values.

Cross-check the vectorized moment formulas against PolyChaos `Statistics.mean` and `Statistics.std` for representative scalar coefficient vectors in tests.

---

## 14. Core-result-specific projection

### 14.1 CableConstants

Fit native:

```text
R
L
C
```

Requirements:

- all collocation results share one concrete CableConstants representation;
- all observed values are finite real numbers;
- coefficient products are vectors with term axis last by definition;
- reconstruct the mean `CableConstants` through its current owner constructor;
- statistics contain native mean and standard deviation for `R`, `L`, and `C`.

### 14.2 LineParameters

Fit native:

```text
R
L
G
C
```

Do not fit:

```text
Z
Y
magnitude
angle
phase
real/imag aliases
```

Requirements:

- all collocation results share the same domain, basis, conductor dimensions, and frequency vector;
- reject zero frequencies because `L` and `C` are undefined through division by angular frequency at zero;
- observed `R/L/G/C` arrays are finite real values;
- coefficient arrays have shape:
  ```text
  conductor × conductor × frequency × basis_term
  ```
- reconstruct the mean core result using:
  ```text
  Z = R + jωL
  Y = G + jωC
  ```
  and the current `LineParameters` owner constructor;
- preserve domain, frequency axis, result basis, ordering, and matrix shape exactly;
- statistics contain native mean and standard deviation arrays for `R`, `L`, `G`, and `C`.

Do not modify the Engine or DataModel core-result representations to accommodate polynomial chaos.

---

## 15. Independent surrogate validation

Validation is mandatory in this implementation.

For each outer point:

1. derive a deterministic local RNG seed from `validation_seed` and the one-based outer-point index without using the global RNG;
2. generate `validation_points` independent latent coordinates from the same supported standard measure:
   - `randn` for normal;
   - uniform `[-1, 1]` for uniform;
3. map them to full physical-value tuples with the same affine maps used for collocation;
4. realise and compute the actual core problem at every validation point through `realize(point, values)`;
5. evaluate the fitted PCE at the latent validation matrix using `PolyChaos.evaluate`/`evaluatePCE` public interfaces and matrix multiplication;
6. compare predicted and actual native `R/L/C[/G]` products separately.

For each quantity family calculate:

```text
relative RMS error = norm(predicted - actual) / max(norm(actual), eps(Float64))
relative max error = maximum(abs(predicted - actual)) /
                     max(maximum(abs(actual)), eps(Float64))
```

The expansion is accepted only when **both** metrics are less than or equal to `validation_rtol` for every fitted quantity.

On failure throw a dedicated owner-defined exception or an `ArgumentError` with all of:

```text
outer point index
distribution
degree
quadrature order
stochastic dimension
basis term count
collocation evaluation count
validation evaluation count
worst quantity
worst matrix/frequency index
relative RMS error
relative maximum error
requested tolerance
```

Do not silently increase degree, refit, retry, or fall back to Monte Carlo.

Store successful validation diagnostics in `validation(result)`.

---

## 16. Observation and publication integration

Integrate with the **live accepted product-aware observation grammar**. Do not create a PolyChaos-only selector system.

At minimum support owner methods equivalent to:

```julia
observe(
    result::PolynomialChaosResult,
    ::typeof(statistics),
    ::typeof(R),
    ::typeof(Statistics.mean),
    point,
    indices...,
)

observe(
    result::PolynomialChaosResult,
    ::typeof(statistics),
    ::typeof(R),
    ::typeof(Statistics.std),
    point,
    indices...,
)
```

and the corresponding `L`, `C`, and `G` methods for the represented core result.

Requirements:

- direct core means remain available through ordinary result-space iteration and existing core-result `observe` methods;
- `observables(result, explicit_requests)` publishes detached values through the existing Units authority;
- PlotBuilder and ReportBuilder are not modified in this campaign beyond compilation fixes required by the established observable protocol;
- no consumer reads `expansion_values`, coefficient arrays, or statistics storage fields directly;
- do not define an envelope-dumping `observables(result)` method;
- do not assign display units during fitting or validation.

No plots or report formats for polynomial-chaos coefficients are added in this campaign.

---

## 17. Supplemental output

When `options.retain_details == false`:

```julia
details(result) == (;)
```

When enabled, retain one fixed-key record per outer `GridPoint`, aligned with result order:

```julia
(
    points = [
        (
            collocation = [...],
            validation = [...],
        ),
        ...
    ],
)
```

Each retained record is produced through the computation owner's existing `computation_details` contract. Dynamic third-party payloads may exist only inside an explicit opaque leaf of that typed record.

Do not place basis objects, coefficients, statistics, validation metrics, counts, or formulation fields into `details`; those are known first-class products.

---

## 18. Future higher-order extension boundary

This campaign must leave the following extension path complete:

```text
new higher-order formulation
    -> dependency-native point generation
    -> points(space)
    -> uncertainties(point)
    -> formulation-owned mapping to concrete physical values
    -> realize(point, values)
    -> compute(core_problem, inner)
    -> formulation-owned accumulation/result
```

The operations divide cleanly:

```text
materialize(point)
    deterministic point construction, including Measurements propagation

realize(rng, point, distribution)
    random stochastic realization used by Monte Carlo

realize(point, values)
    supplied-value realization used by deterministic stochastic designs
    such as polynomial-chaos collocation, QMC points, or sparse grids
```

Adding future quasi-Monte Carlo, Latin hypercube, sparse-grid, sigma-point, or related methods must not require changing:

- `Grid`;
- `GridPoint`;
- `Gridspace`;
- `UncertainValue`;
- CableDesign/System builders;
- Engine kernels;
- `uncertainties`;
- the existing meanings of `materialize` or `realize`;
- a central sampler registry;
- a realization-policy hierarchy.

A future formulation obtains coordinates through the dependency that owns them, maps those coordinates to a concrete physical-value tuple, and calls the existing `realize(point, values)` method. Do not add QuasiMonteCarlo, Sobol sampling, sparse grids, adaptive PCE, regression PCE, or sensitivity indices in this campaign.

---

## 19. Tests

Use existing test infrastructure and extend current files by owner. Do not create a parallel test framework.

### 19.1 ParametricBuilder contracts

Test:

- existing deterministic materialisation remains unchanged;
- existing Measurements-backed `UncertainValue` materialisation remains unchanged;
- existing RNG-backed stochastic realisation remains unchanged;
- deterministic uncertainty ordering for flat and nested `GridPoint`s;
- one builder argument reused internally occupies one coordinate;
- equal descriptors in two slots occupy two coordinates;
- zero-sigma descriptors remain in the positional tuple;
- `realize(point, values)` reconstructs the expected target;
- nested callables are invoked child-first as before;
- cardinality mismatch fails before builder invocation;
- invalid supplied physical values surface owner validation;
- representative calls are inferred;
- scalar-target supplied-value realisation adds zero allocations after compilation;
- `realize` remains non-root-exported but available to owner extensions;
- no `Realization*` policy/type hierarchy exists;
- Monte Carlo fixed-seed outputs, draw order, seeds, trial counts, statistics, retained products, histograms, and details remain unchanged.

### 19.2 Extension unloaded

Without PolyChaos loaded:

- LineCableModels loads and precompiles;
- the PolyChaos extension is absent;
- formulation construction succeeds;
- computation has no fake fallback.

### 19.3 PolyChaos dependency API

Against the resolved dependency version, test:

- requested univariate quadrature point count;
- `MultiOrthoPoly` basis term count `binomial(d + p, p)`;
- tensor-node count `q^d`;
- tensor node ordering and product weights;
- basis-evaluation matrix orientation;
- affine mean/std mapping for normal and uniform descriptors;
- vectorized coefficient projection against a scalar reference calculation;
- vectorized mean/std against PolyChaos's public statistics functions.

### 19.4 Exact analytic core problems

Create real test-only LineCableModels problem/formulation methods, not mocks.

The test solver returns `CableConstants` or `LineParameters` whose native `R/L/C/G` values are known polynomials of two uncertain inputs.

Required cases:

- degree-1 normal response reproduced to numerical tolerance;
- degree-2 normal response reproduced to numerical tolerance;
- degree-1 uniform response reproduced to numerical tolerance;
- analytic mean and standard deviation match the fitted results;
- independent validation succeeds;
- an intentionally under-degree nonlinear response fails validation;
- validation error reports the required diagnostics.

### 19.5 Real LineCableModels smoke test

Use one small actual cable/line-parameter problem with two or three independent uncertain material properties and a low degree/order that stays cheap in CI.

Prove:

- result core mean type and metadata;
- coefficient and statistic shapes;
- finite validation metrics;
- fewer actual core solves than a modest Monte Carlo reference configured for the test;
- no change to deterministic core results.

Do not turn a slow research benchmark into a mandatory CI test.

### 19.6 Cost and failure boundaries

Test:

- total cost rejection occurs before the first core solve;
- unsupported distribution fails at construction;
- correlated-input requests are not silently accepted;
- zero active dimension fails before solves;
- invalid physical quadrature values surface the owner constructor/validation error;
- collocation result type, frequency, basis, domain, and shape mismatches fail explicitly;
- zero frequency in LineParameters fails explicitly;
- nonfinite core observables fail explicitly;
- unsupported core result family raises `MethodError` or the repository-standard unsupported-operation error without returning an empty result.

### 19.7 Result and observables

Test:

- `PolynomialChaosResult` collection interface and `eltype`;
- `only` works through Base for singleton result spaces;
- `statistics`, `expansions`, `validation`, and `details` alignment;
- explicit observation/publication of mean and standard deviation for every supported quantity;
- Units conversion happens only at publication;
- missing unsupported products such as `samples` and `histograms` are not fabricated;
- method ambiguity and inference checks pass.

---

## 20. Benchmark and documentation

Add one executable, non-CI benchmark under the repository's established benchmark or development owner, comparing:

```text
Monte Carlo reference
PolynomialChaos degree 1
PolynomialChaos degree 2
PolynomialChaos degree 3
PolynomialChaos degree 4
```

Use a small cable case with three uncertain material properties. Record:

```text
actual core solve count
wall time
mean R/L/G/C error
standard-deviation error
independent validation error
selected surrogate quantiles computed offline from the fitted expansion
```

Do not expose surrogate samples through `samples(result)`.

Update maintained documentation with:

- dependency activation through `using PolyChaos`;
- the exact public constructor/call;
- normal and uniform semantics;
- independence restriction;
- tensor cost law `q^d`;
- hard evaluation budget;
- mandatory validation;
- why `R/L/C/G` are fitted instead of `Z/Y` magnitude/angle;
- distinction between core-solver trials and cheap surrogate evaluation;
- one complete example for CableConstants or LineParameters;
- future extension contract using `uncertainties`, `materialize`, and `realize`.

Do not add marketing claims that PCE is universally faster. State that it is appropriate for low-dimensional smooth responses and is rejected by validation or budget guards otherwise.

---

## 21. Commit sequence

Execute these five commits. Keep every commit green and reviewable.

### Commit 1

Subject:

```text
feat(parametricbuilder): extend uncertainty realisation
```

Scope:

- preserve and document the existing materialisation/realisation distinction;
- add `uncertainties(point)`;
- add `realize(point, values)`;
- make the internal owner-visible ParametricBuilder seam explicit;
- retain existing `materialize` and RNG-backed `realize` methods unchanged unless a behavior-identical local refactor is proven by fixed-seed and allocation tests;
- add focused ordering, nesting, cardinality, inference, allocation, Measurements, and Monte Carlo preservation tests.

Forbidden:

- PolyChaos dependency;
- PolynomialChaos formulation/result types;
- removal, renaming, aliasing, or root export of `realize`;
- realization-policy types or registries;
- Monte Carlo numerical or API changes.

### Commit 2

Subject:

```text
feat(uq): define polynomial chaos contracts
```

Scope:

- add `PolynomialChaos`;
- add `PolynomialChaosResult`;
- add `expansions` and `validation` accessors;
- add constructor/result invariants and core tests;
- update package/module exports.

Forbidden:

- importing PolyChaos in core;
- fake compute fallback;
- plotting/reporting features.

### Commit 3

Subject:

```text
feat(polychaos): add spectral projection extension
```

Scope:

- add weak dependency, compat, extension, extras, and target;
- implement basis construction, tensor quadrature, physical mapping, cost preflight, collocation, vectorized fitting, mean reconstruction, statistics, validation, and details;
- add focused extension tests in the same commit so it remains green.

Forbidden:

- adaptive degree;
- sparse grids;
- regression PCE;
- QuasiMonteCarlo;
- core solver changes.

### Commit 4

Subject:

```text
test(polychaos): validate accuracy and integration
```

Scope:

- add exact analytic CableConstants and LineParameters tests;
- add real solver smoke test;
- add budget/failure/inference/ambiguity tests;
- add result-observation integration tests.

Forbidden:

- production API redesign to make tests easier;
- tolerance relaxation without numerical justification.

### Commit 5

Subject:

```text
docs(uq): document polynomial chaos workflow
```

Scope:

- add user/developer documentation;
- add the non-CI benchmark;
- update extension and API references;
- perform final residue and dependency-boundary audit.

Forbidden:

- unrelated cleanup;
- adding another UQ method.

---

## 22. Validation

Discover and run the live repository's exact commands. At minimum run the applicable equivalents of:

```bash
julia --project=. --color=yes -e 'using Pkg; Pkg.test()'
git diff --check
```

Also run:

- focused ParametricBuilder tests;
- focused UQ tests;
- core-only extension-unloaded tests;
- PolyChaos extension tests;
- exact analytic PCE tests;
- real CableConstants/LineParameters smoke tests;
- Aqua;
- method ambiguity checks;
- formatter/linter gates;
- documentation doctests and build;
- package precompilation both with and without PolyChaos loaded;
- scoped residue searches proving that `realize` remains, retired alternative spellings are absent, and no realization-policy hierarchy was introduced.

Inspect representative concrete paths with:

```text
@inferred
@code_warntype
@allocated
```

Do not weaken existing numerical tolerances, regenerate unrelated baselines, skip failing tests, or silently exclude the PolyChaos extension from the full suite.

---

## 23. Hard prohibitions

Do not:

- provide another plan instead of implementing;
- remove, rename, bypass, deprecate, or compatibility-wrap the existing `realize` generic;
- introduce a second realisation verb under another spelling;
- add a universal UQ sampler/generator/point-design abstraction;
- add a realization policy hierarchy;
- add identity keys, caches, configurations, bindings, or manifests to Gridspace;
- change Gridspace product/zip semantics;
- put PolyChaos in core `[deps]`;
- import PolyChaos from core source;
- make the deterministic Engine aware of polynomial chaos;
- implement intrusive stochastic Galerkin equations;
- fit display-unit values;
- fit complex magnitude or wrapped phase;
- clamp invalid physical quadrature values;
- silently normalize malformed dependency output;
- use PolyChaos private APIs;
- create dense diagonal `H` or `W` matrices;
- perform one scalar regression per matrix/frequency entry;
- auto-increase degree;
- silently fall back to Monte Carlo;
- exceed `max_evaluations` after warning;
- fabricate empirical samples or histograms from PCE and expose them as solver trials;
- add covariance support without an explicit independent-latent transform;
- add QMC, Sobol sensitivity, sparse grids, adaptive PCE, or surrogate plotting;
- introduce `helpers.jl`, `utils.jl`, `common.jl`, managers, registries, or pass-through wrappers;
- import or extend another owner's underscore-prefixed method;
- touch unrelated PlotBuilder, ReportBuilder, Units, XLSX, or visualization behavior.

---

## 24. Completion criteria

The goal is complete only when all of the following hold:

1. Existing `materialize(point)` semantics remain unchanged.
2. Existing `realize(rng, point, distribution)` semantics and fixed-seed behavior remain unchanged.
3. `uncertainties(point)` defines one deterministic structural coordinate order matching recursive realisation.
4. `realize(point, values)` constructs the complete target from supplied concrete physical values in that order.
5. Materialisation and realisation remain internal, non-root-exported ParametricBuilder protocols.
6. No competing realisation verb, mapper, policy, context, registry, or hierarchy exists.
7. PolyChaos is a weak dependency activated through a Julia package extension.
8. Core LineCableModels loads and precompiles without PolyChaos.
9. `PolynomialChaos` is a first-class UQ formulation parallel to MonteCarlo and LinearError.
10. `PolynomialChaosResult` is a proper result space aligned with outer Gridspace order.
11. Only CableConstants and LineParameters are supported initially.
12. Normal and uniform independent uncertainty semantics match the declared mean and standard uncertainty.
13. Tensor quadrature size and total core-solve cost are known and checked before execution.
14. Collocation nodes are mapped to full physical-value tuples and passed through `realize(point, values)`.
15. Spectral coefficients are fitted vectorially in native R/L/C[/G] values.
16. Mean core results preserve domain, basis, frequency, shape, and ordering.
17. Statistics derive from coefficients rather than surrogate resampling.
18. Independent validation is mandatory and failures are explicit.
19. `statistics`, `expansions`, `validation`, `details`, iteration, `eltype`, and Base `only` work as specified.
20. Explicit observation/publication of PCE means and standard deviations uses the existing Grammar and Units authorities.
21. No fake `samples` or `histograms` product exists.
22. Exact polynomial tests, real solver smoke tests, cost guards, extension boundaries, inference, ambiguity, docs, and full tests pass.
23. No new generic sampler bureaucracy, Gridspace administration, core-solver stochastic code, or unrelated refactor has been introduced.

---

## 25. Final execution report

After implementation, report only:

- starting branch and HEAD;
- five commit hashes and subjects;
- production files changed per commit;
- tests/docs/benchmark files changed per commit;
- exact PolyChaos version resolved;
- exact validation commands and outcomes;
- core-solve counts for analytic and real smoke cases;
- confirmation that fixed-seed Monte Carlo behavior was preserved;
- confirmation that existing RNG-backed `realize` behavior is unchanged and the supplied-value overload is covered;
- confirmation that no competing realisation vocabulary or policy hierarchy was introduced;
- confirmation that LineCableModels loads with PolyChaos absent and extension loads with PolyChaos present;
- final public API example;
- any stop-condition evidence that prevented completion.

Do not append speculative follow-up architecture.