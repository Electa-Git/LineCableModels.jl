_validate_shunt_policy(::AbstractFormulation) = nothing
function _validate_shunt_policy(inner::Union{Engine.LineParametersFormulation,Engine.CableConstantsFormulation})
    selected = inner.methods.shunt_model
    if selected isa Engine.ShuntModel.Formula{:boundary} && selected.parameters.fallback !== :error
        throw(ArgumentError("uncertainty propagation requires a fixed shunt model; select strict :boundary or :coaxial, without automatic fallback"))
    end
    return nothing
end

"""
$(TYPEDEF)

Select direct linear uncertainty propagation with `inner`.

Shared inputs must be supplied once to a joint `Gridspace` builder, which
derives the dependent geometry. Propagation retains those correlations and
differentiates continuous geometry at the nominal design, including bounded
power-cell compaction. Strand counts and clipping topology are nominal discrete
choices: a study crossing such a transition has no single smooth linear model.
First-order moments need not equal nonlinear Monte Carlo moments.

$(TYPEDFIELDS)
"""
struct LinearError{F <: AbstractFormulation, O <: ComputationOptions} <: AbstractFormulation
    "Formulation used for each materialized problem."
    inner::F
    "Supplemental-output retention options owned by this propagation."
    options::O

    function LinearError(inner::AbstractFormulation, options::ComputationOptions)
        _validate_shunt_policy(inner)
        normalized = computation_options(LinearError, options)
        return new{typeof(inner), typeof(normalized)}(inner, normalized)
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

function LinearError(
        inner::F;
        options::Union{NamedTuple,ComputationOptions} = ComputationOptions()
) where {F <: AbstractFormulation}
    options = options isa NamedTuple ? ComputationOptions(options) : options
    return LinearError(inner, options)
end

function Base.NamedTuple(value::LinearError)
    return (kind = :linear_error, inner = NamedTuple(value.inner), options = value.options.data)
end

"""Identify first-order uncertainty propagation without inspecting its inner owner."""
description(::Type{<:LinearError}; compact::Bool=false) = "LEP"
description(::LinearError; compact::Bool=false) = description(LinearError;compact)
formula_id(::Type{<:LinearError}) = :LinearError
formula_id(::LinearError) = :LinearError
computation_options(value::LinearError) = value.options
function Base.pairs(value::LinearError; quantity = nothing)
    pairs(LinearError, (inner = value.inner, options = value.options.data); quantity)
end

"""
$(TYPEDEF)

Select conditional Monte Carlo propagation over a
[`ParametricProblem`](@ref). Randomness is local and reproducible when `seed`
is supplied. Computation option `on_error=:fail` propagates every exception.
`on_error=:retry` rejects only realizations that raise `DomainError`, retains
their sampled arguments and error summaries, and continues until the requested
number of successful trials is obtained or `max_failures` is reached. Retry
mode requires `retain_details=true` and estimates the output distribution
conditional on successful problem construction and computation.

Load `Measurements` before computation. Aggregation stores marginal
uncertainty-bearing cores from accepted sample means and sample standard
deviations, not histogram bins or standard errors of the mean. Sampling and
histogram retention do not change this payload. Joint output correlations are
not retained by the marginal surrogate. Access and transport reuse the stored
uncertainty-source identities.

Line-system construction enforces exterior clearance on every realization.
When Measurements is loaded, the propagated clearance reserve is prepared
before sampling and retained for every draw. Adjusted placements therefore
describe a clearance-constrained system. Adjustments emit one summary per
parameter point, not one warning per trial; retained details include their
count and maximum displacement \\[m\\]. Invalid dimensions and infeasible
internal cable constructions remain errors subject to `on_error`.

Use one joint `Gridspace` builder for dependent internal dimensions; derive
stack boundaries and wire placements from its shared inputs. Repeating a Grid
as separate source arguments creates independent draws. The default normal
law has unbounded support: marginal means and standard uncertainties cannot
guarantee feasible geometry. `distribution=:uniform` instead has support
`nominal ± sqrt(3)*sigma` for each primitive descriptor; the builder must map
the complete joint support to feasible designs. Selecting this law is an
explicit statistical assumption, not a repair or an inferred correlation.

$(TYPEDFIELDS)
"""
struct MonteCarlo{F <: AbstractFormulation, O <: ComputationOptions} <:
       AbstractFormulation
    "Formulation used for each sampled problem."
    inner::F
    "Normalized sampling, error-handling and supplemental-output computation options."
    options::O

    function MonteCarlo(inner::AbstractFormulation, options::ComputationOptions)
        _validate_shunt_policy(inner)
        normalized = computation_options(MonteCarlo, options)
        return new{typeof(inner), typeof(normalized)}(inner, normalized)
    end
end

"""Identify Monte Carlo propagation without sampling or configuring a solver."""
description(::Type{<:MonteCarlo}; compact::Bool=false) = "Monte Carlo"
description(::MonteCarlo; compact::Bool=false) = description(MonteCarlo;compact)
description(::Type{<:LinearError},::Val{:representation};compact::Bool=false) =
    "dependency-preserving uncertainty"
description(::Type{<:MonteCarlo},::Val{:representation};compact::Bool=false) =
    "marginal mean ± std"

# Capture display annotations separately from the estimator/representation IDs
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
    if ismissing(retained.inner)
        push!(entries, (owner, (:inner,)) => missing)
    else
        append!(entries,
            pairs((retained.inner isa Pair ? Tuple(retained.inner) : (retained.inner,))...; quantity))
    end
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
    normalized.on_error === :fail || normalized.on_error === :retry ||
        throw(ArgumentError(
            "MonteCarlo on_error must be :fail or :retry",
        ))
    normalized.on_error === :retry && !normalized.retain_details &&
        throw(
            ArgumentError(
            "MonteCarlo on_error=:retry requires retain_details=true",
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

Construct a Monte Carlo calculation with execution controls supplied as an
ordinary `options` named tuple or as keyword shorthand. Both forms are
normalized by `computation_options(MonteCarlo, options)`. A key supplied in
both places is an error.

# Arguments

- `inner`: Formulation used for each sampled problem.

# Keywords

`options=(;)` holds any of the following controls. Additional keywords use
the same names and are merged into `options` before normalization.

- `trials=nothing`: Positive accepted-trial count, or DKW sizing when omitted.
- `confidence=0.95`: Simultaneous empirical-CDF confidence [dimensionless], in `(0, 1)`.
- `cdf_tol=0.02`: Maximum empirical-CDF deviation for DKW sizing [dimensionless], in `(0, 1)`.
- `distribution=:normal`: `:normal`, `:uniform`, a sampler function, or an extension-supported distribution.
- `seed=nothing`: Nonnegative root seed representable as `UInt64`, or fresh randomness.
- `return_samples=false`: Retain joint samples.
- `return_histograms=false`: Retain marginal histogram densities.
- `bins=nothing`: Positive histogram bin count, or automatic binning.
- `retain_details=false`: Retain accepted-trial details and failure diagnostics.
- `on_error=:fail`: Propagate exceptions; `:retry` rejects `DomainError` realizations and requires `retain_details=true`.
- `max_failures=100`: Positive maximum rejected-trial count per parameter point.

# Returns

- A `MonteCarlo` calculation storing normalized `ComputationOptions` in `options`.
  The concrete payload type retains the sampler type.
"""
function MonteCarlo(inner::AbstractFormulation; options::Union{NamedTuple,ComputationOptions} = ComputationOptions(), kwargs...)
    options = options isa NamedTuple ? ComputationOptions(options) : options
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
