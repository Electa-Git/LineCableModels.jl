"""
$(TYPEDEF)

Select the derivative calculation used by [`LinearError`](@ref).
"""
abstract type LinearErrorMethod end

"""
$(TYPEDEF)

Propagate first-order input dependencies through Measurements arithmetic in one
computation. Numerical operations must support Measurements scalars.
"""
struct MeasurementsPropagation <: LinearErrorMethod end

"""
$(TYPEDEF)

Estimate the Jacobian by paired central perturbations through the Gridspace
builder. Each perturbation is `relative_step` times the input standard
uncertainty. The default relative step is `0.01`.

$(TYPEDFIELDS)
"""
struct CentralDifferences <: LinearErrorMethod
    "Ratio of the perturbation to the input standard uncertainty [dimensionless]."
    relative_step::Float64

    function CentralDifferences(relative_step)
        validate(relative_step, CentralDifferences)
        return new(Float64(relative_step))
    end
end

CentralDifferences(; relative_step=0.01) = CentralDifferences(relative_step)

# Public symbols are resolved here, before the computation is constructed.
function LinearErrorMethod(identifier::Symbol, relative_step=0.01)
    if identifier === :measurements
        validate(relative_step, CentralDifferences)
        return MeasurementsPropagation()
    elseif identifier === :central
        return CentralDifferences(relative_step)
    end
    throw(ArgumentError("LinearError derivative must be :measurements or :central"))
end
LinearErrorMethod(method::LinearErrorMethod) = method
function LinearErrorMethod(method::LinearErrorMethod, relative_step)
    throw(ArgumentError("relative_step=$(repr(relative_step)) was supplied alongside " *
        "$(typeof(method)); specify the step in CentralDifferences"))
end
function LinearErrorMethod(value, args...)
    throw(ArgumentError("LinearError derivative must be :measurements, :central or " *
        "a LinearErrorMethod value; got $(repr(value))"))
end

function validate(step::Real, ::Type{CentralDifferences})
    isfinite(step) && step > 0 && isfinite(Float64(step)) && Float64(step) > 0 ||
        throw(ArgumentError("LinearError relative_step must be finite and positive"))
    return step
end
function validate(step::Bool, ::Type{CentralDifferences})
    throw(ArgumentError("LinearError relative_step must be a real number, excluding Bool; got $step"))
end
function validate(step, ::Type{CentralDifferences})
    throw(ArgumentError("LinearError relative_step must be a finite positive real number; got $(repr(step))"))
end

# Symbols here describe the method for display and serialization only.
Base.NamedTuple(::MeasurementsPropagation) = (derivative=:measurements, relative_step=0.01)
Base.NamedTuple(method::CentralDifferences) = (derivative=:central, relative_step=method.relative_step)

"""
$(TYPEDEF)

Select first-order uncertainty propagation with `inner`.

Option `derivative=:measurements` (default) propagates Measurements through one
evaluation and requires every numerical operation to support those scalars.
Option `derivative=:central` evaluates the nominal problem and two perturbations
per independent uncertain input through the same Gridspace builder. It supports
modal and transient computations with ordinary numerical scalars. Both methods
need `using Measurements` to construct uncertainty-bearing results.

The symbols select UQ-owned method values at construction. The equivalent direct
forms are `derivative=MeasurementsPropagation()` and
`derivative=CentralDifferences(relative_step=0.01)`, imported from
`LineCableModels.UQ`. A method value stores its configuration. Set the relative
step in the `CentralDifferences` constructor.

For central differences, `relative_step=0.01` sets the dimensionless ratio
``s = h_i/u_i`` of each perturbation to its input standard uncertainty. The
Jacobian and output standard uncertainty are

```math
J_i = \\frac{y(x+h_i e_i)-y(x-h_i e_i)}{2h_i}, \\qquad
u_y = \\sqrt{\\sum_i (J_i u_i)^2}.
```

Inputs declared by Grid are independent. Shared dimensions derived by one joint
builder keep their dependence. Central differences reject supplied uncertain
scalars whose shared or correlated sources are outside that declaration.

Modal perturbations use nominal mode assignment and report overlaps in
`details(result).data.points`. A rejected perturbation throws an error whose
`details(error)` identifies the input path, signed step and original failure.
The method preserves the requested step. `retain_details=true` also keeps the scalar
computation records for all evaluations.

Shared inputs must be supplied once to a joint `Gridspace` builder, which
derives the dependent geometry. Propagation retains those correlations and
differentiates continuous geometry at the nominal design, including bounded
power-cell compaction. Strand counts and clipping topology are nominal discrete
choices. A smooth linear model cannot describe a study across such a transition.
First-order moments need not equal nonlinear Monte Carlo moments.

$(TYPEDFIELDS)
"""
struct LinearError{F <: AbstractFormulation, O <: ComputationOptions,
        D <: LinearErrorMethod} <: AbstractFormulation
    "Formulation used for each materialized problem."
    inner::F
    "Derivative method and its numerical configuration."
    derivative::D
    "Supplemental-output retention options owned by this propagation."
    options::O

    function LinearError(inner::AbstractFormulation, derivative::LinearErrorMethod,
            options::ComputationOptions)
        validate(inner, LinearError)
        normalized = computation_options(LinearError, options)
        return new{typeof(inner), typeof(normalized), typeof(derivative)}(
            inner, derivative, normalized)
    end
end

function computation_options(
        ::Type{LinearError},
        record::ComputationOptions
)::ComputationOptions
    options = record.data
    unknown = filter(key -> key !== :retain_details, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown LinearError computation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge((retain_details = false,), options)
    normalized.retain_details isa Bool || throw(ArgumentError(
        "LinearError retain_details must be Bool",
    ))
    return ComputationOptions(; retain_details = normalized.retain_details)
end

function LinearError(inner::AbstractFormulation, options::ComputationOptions)
    selected = get(options.data, :derivative, :measurements)
    derivative = haskey(options.data, :relative_step) ?
        LinearErrorMethod(selected, options.data.relative_step) : LinearErrorMethod(selected)
    unknown = filter(key -> key ∉ (:retain_details,:derivative,:relative_step), keys(options.data))
    isempty(unknown) || throw(ArgumentError(
        "unknown LinearError computation options: $(sort!(collect(unknown)))"))
    retention = ComputationOptions(; retain_details=get(options.data,:retain_details,false))
    return LinearError(inner, derivative, retention)
end

function LinearError(inner::AbstractFormulation, derivative::LinearErrorMethod;
        options::Union{NamedTuple,ComputationOptions}=ComputationOptions())
    return LinearError(inner, derivative, ComputationOptions(options))
end

function LinearError(
        inner::F;
        options::Union{NamedTuple,ComputationOptions} = ComputationOptions(),
        kwargs...
) where {F <: AbstractFormulation}
    options = ComputationOptions(options)
    supplied = (; kwargs...)
    duplicates = intersect(keys(options.data), keys(supplied))
    isempty(duplicates) || throw(ArgumentError(
        "LinearError options were supplied twice: $(collect(duplicates))"))
    return LinearError(inner, ComputationOptions(merge(options.data, supplied)))
end

function Base.NamedTuple(value::LinearError)
    return (kind = :linear_error, inner = NamedTuple(value.inner), options = merge(value.options.data, NamedTuple(value.derivative)))
end

"""Identify first-order uncertainty propagation without inspecting its inner owner."""
description(::Type{<:LinearError}; compact::Bool=false) = "LEP"
description(::LinearError; compact::Bool=false) = description(LinearError;compact)
formula_id(::Type{<:LinearError}) = :LinearError
formula_id(::LinearError) = :LinearError
computation_options(value::LinearError) = value.options
function Base.pairs(value::LinearError; quantity = nothing)
    pairs(LinearError, (inner = value.inner,
        options = merge(value.options.data, NamedTuple(value.derivative))); quantity)
end

"""
$(TYPEDEF)

Select conditional Monte Carlo propagation over a
[`ParametricProblem`](@ref). Randomness is local and reproducible when `seed`
is supplied. Computation option `on_error=:fail` propagates every exception.
`on_error=:resample` rejects only realizations that raise `DomainError` and retains
their sampled arguments and error summaries. Sampling continues until the requested
number of successful trials is obtained or `max_failures` is reached. Resampling
mode requires `retain_details=true` and estimates the output distribution
conditional on successful problem construction and computation.

Load `Measurements` before computation. Aggregation stores marginal
uncertainty-bearing cores from accepted sample means and sample standard
deviations, not histogram bins or standard errors of the mean. This payload is independent of sampling and histogram retention. Joint output correlations are
not retained by the marginal surrogate. Access and transport reuse the stored
uncertainty-source identities.

Line-system construction checks exterior clearance on every realization. A draw
whose cables overlap throws
[`OverlapError`](@ref LineCableModels.DataModel.OverlapError) at its build. The
study discards that draw as if it was never drawn and draws again, whatever
`on_error` says. It records nothing about the discarded draw. `max_failures`
overlapping draws in a row stop the study. Invalid dimensions and infeasible internal
cable constructions remain errors subject to `on_error`.

Use one joint `Gridspace` builder for dependent internal dimensions. Derive
stack geometric boundaries and wire placements from its shared inputs. Repeating a Grid
as separate source arguments creates independent draws. The default normal
law has unbounded support: marginal means and standard uncertainties cannot
guarantee feasible geometry. `distribution=:uniform` instead has support
`nominal ± sqrt(3)*sigma` for each primitive descriptor. The builder must map
the complete joint support to feasible designs. Select this law only when the
uniform distribution is a justified statistical assumption. Changing the
marginal law neither repairs the geometry nor supplies a correlation.

$(TYPEDFIELDS)
"""
struct MonteCarlo{F <: AbstractFormulation, O <: ComputationOptions} <:
       AbstractFormulation
    "Formulation used for each sampled problem."
    inner::F
    "Normalized sampling, error-handling and supplemental-output computation options."
    options::O

    function MonteCarlo(inner::AbstractFormulation, options::ComputationOptions)
        validate(inner, MonteCarlo)
        normalized = computation_options(MonteCarlo, options)
        return new{typeof(inner), typeof(normalized)}(inner, normalized)
    end
end

# Uncertainty propagation accepts an inner formulation whose shunt model is fixed.
validate(inner::AbstractFormulation, ::Type{<:Union{LinearError, MonteCarlo}}) = inner
function validate(inner::Union{Engine.LineParametersFormulation, Engine.CableConstantsFormulation},
        propagation::Type{<:Union{LinearError, MonteCarlo}})
    validate(inner.formulas.shunt_model, propagation)
    return inner
end
function validate(selected::Engine.ShuntModel.Formula{:boundary},
        ::Type{<:Union{LinearError, MonteCarlo}})
    selected.parameters.fallback === :error || throw(ArgumentError("uncertainty propagation requires a fixed shunt model; select strict :boundary or :equivalent, without automatic fallback"))
    return selected
end

"""Identify Monte Carlo propagation without sampling or configuring a solver."""
description(::Type{<:MonteCarlo}; compact::Bool=false) = "Monte Carlo"
description(::MonteCarlo; compact::Bool=false) = description(MonteCarlo;compact)
description(::Type{<:LinearError},::Val{:representation};compact::Bool=false) =
    "dependency-preserving uncertainty"
description(::Type{<:MonteCarlo},::Val{:representation};compact::Bool=false) =
    "marginal mean ± std"

# Capture display annotations separately from the estimator and representation IDs
# used by scientific grouping. Both names remain supplied by their UQ owner.
function _uncertainty_descriptions(owner)
    return (estimator=(name="uncertainty estimator",unit="",text=description(owner;compact=true)),
        representation=(name="uncertainty representation",unit="",
            text=description(owner,Val(:representation);compact=true)))
end
formula_id(::Type{<:MonteCarlo}) = :MonteCarlo
formula_id(::MonteCarlo) = :MonteCarlo
computation_options(value::MonteCarlo) = value.options
function Base.pairs(value::MonteCarlo; quantity = nothing)
    pairs(MonteCarlo, (inner = value.inner, options = value.options.data); quantity)
end
function Base.pairs(owner::Type{<:Union{MonteCarlo, LinearError}}, retained::NamedTuple; quantity = nothing)
    entries=Pair{Tuple, Any}[(owner, ()) => (owner => retained.options)]
    append!(entries, pairs(retained.inner; quantity))
    return entries
end
description(::Type{<:Union{MonteCarlo, LinearError}}, ::Val{:inner}; compact::Bool=false) = "inner method"

function computation_options(
        ::Type{MonteCarlo},
        record::ComputationOptions
)::ComputationOptions
    options = record.data
    defaults = (
        trials = nothing, confidence = 0.95, cdf_tol = 0.02,
        distribution = :normal, seed = nothing,
        return_samples = false, return_histograms = false, bins = nothing,
        retain_details = false, on_error = :fail, max_failures = 100
    )
    unknown = filter(key -> key ∉ keys(defaults), keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown MonteCarlo computation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge(defaults, options)
    for key in (:trials, :bins, :max_failures)
        value = getproperty(normalized, key)
        key !== :max_failures && value === nothing && continue
        value isa Integer && !(value isa Bool) && 0 < value <= typemax(Int) ||
            throw(ArgumentError("MonteCarlo $key must be a positive integer representable as Int"))
    end
    for key in (:confidence, :cdf_tol)
        value = getproperty(normalized, key)
        value isa Real && !(value isa Bool) && 0 < value < 1 &&
        0 < Float64(value) < 1 || throw(ArgumentError(
            "MonteCarlo $key must lie strictly between zero and one as Float64",
        ))
    end
    for key in (:return_samples, :return_histograms, :retain_details)
        getproperty(normalized, key) isa Bool || throw(ArgumentError(
            "MonteCarlo $key must be Bool",
        ))
    end
    seed = normalized.seed
    seed === nothing ||
        (seed isa Integer && !(seed isa Bool) && 0 <= seed <= typemax(UInt64)) ||
        throw(ArgumentError("MonteCarlo seed must be a nonnegative integer representable as UInt64"))
    distribution = normalized.distribution
    distribution isa Symbol && distribution ∉ (:normal, :uniform) &&
        throw(ArgumentError(
            "unsupported distribution $(repr(distribution)); expected :normal, :uniform, a sampler function, or an extension-supported distribution",
        ))
    normalized.on_error === :fail || normalized.on_error === :resample ||
        throw(ArgumentError(
            "MonteCarlo on_error must be :fail or :resample",
        ))
    normalized.on_error === :resample && !normalized.retain_details &&
        throw(
            ArgumentError(
            "MonteCarlo on_error=:resample requires retain_details=true",
        ),
        )
    return ComputationOptions(;
        trials = normalized.trials === nothing ? nothing : Int(normalized.trials),
        confidence = Float64(normalized.confidence),
        cdf_tol = Float64(normalized.cdf_tol),
        distribution = distribution,
        seed = seed === nothing ? nothing : UInt64(seed),
        return_samples = normalized.return_samples,
        return_histograms = normalized.return_histograms,
        bins = normalized.bins === nothing ? nothing : Int(normalized.bins),
        retain_details = normalized.retain_details,
        on_error = normalized.on_error,
        max_failures = Int(normalized.max_failures)
    )
end

"""
$(TYPEDSIGNATURES)

Construct a Monte Carlo formulation with execution controls supplied as an
ordinary `options` named tuple or as keyword shorthand. Both forms are
normalized by `computation_options(MonteCarlo, options)`. A key supplied in
both places is an error.

# Arguments

- `inner`: formulation used for each sampled problem.

# Keywords

`options=(;)` accepts these controls. Additional keywords use
the same names and are merged into `options` before normalization.

- `trials=nothing`: positive accepted-trial count, or DKW sizing when omitted.
- `confidence=0.95`: simultaneous empirical-CDF confidence [dimensionless], in `(0, 1)`.
- `cdf_tol=0.02`: maximum empirical-CDF deviation for DKW sizing [dimensionless], in `(0, 1)`.
- `distribution=:normal`: `:normal`, `:uniform`, a sampler function, or an extension-supported distribution.
- `seed=nothing`: nonnegative root seed representable as `UInt64`, or fresh randomness.
- `return_samples=false`: retain joint samples.
- `return_histograms=false`: retain marginal histogram densities.
- `bins=nothing`: positive histogram bin count, or automatic binning.
- `retain_details=false`: retain accepted-trial details and failure diagnostics.
- `on_error=:fail`: propagate exceptions. `:resample` rejects `DomainError` realizations and requires `retain_details=true`.
- `max_failures=100`: positive maximum rejected-trial count per parameter point.

# Returns

- A `MonteCarlo` formulation storing normalized `ComputationOptions` in `options`.
  The concrete payload type retains the sampler type.
"""
function MonteCarlo(inner::AbstractFormulation; options::Union{NamedTuple,ComputationOptions} = ComputationOptions(), kwargs...)
    options = ComputationOptions(options)
    supplied = (; kwargs...)
    duplicates = filter(key -> haskey(options.data, key), keys(supplied))
    isempty(duplicates) || throw(ArgumentError(
        "MonteCarlo computation options supplied both as keywords and in options: $(collect(duplicates))",
    ))
    return MonteCarlo(inner, isempty(supplied) ? options : ComputationOptions(merge(options.data, supplied)))
end

function Base.NamedTuple(value::MonteCarlo)
    return (kind = :monte_carlo, inner = NamedTuple(value.inner), options = value.options.data)
end
