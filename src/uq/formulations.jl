"""
$(TYPEDEF)

Select direct linear uncertainty propagation with `inner`.

Shared inputs must be supplied once to a joint `Gridspace` builder, which
derives the dependent geometry. Propagation retains those correlations and
differentiates continuous geometry at the nominal design, including bounded
power-cell compaction. Strand counts and clipping topology are nominal discrete
choices. A smooth linear model cannot describe a study across such a transition.
First-order moments need not equal nonlinear Monte Carlo moments.

$(TYPEDFIELDS)
"""
struct LinearError{F <: AbstractFormulation, O <: ComputationOptions} <: AbstractTraversal
    "Formulation used for each materialized problem."
    inner::F
    "Supplemental-output retention options owned by this propagation."
    options::O

    LinearError(inner::F, options::O) where {F <: AbstractFormulation, O <: ComputationOptions} =
        validate(new{F, O}(inner, options))
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
is supplied.

A draw is infeasible when its build throws
[`OverlapError`](@ref LineCableModels.DataModel.OverlapError), or when its build or
computation throws `DomainError`. The study discards an infeasible draw as if it was never
drawn and draws again. It records nothing about the discarded draw. After `max_failures`
infeasible draws in a row for one trial, the study throws an error with that count and the
last exception. Any other exception propagates. The output distribution is conditional on
feasible draws.

Load `Measurements` before computation. Aggregation stores marginal
uncertainty-bearing cores from accepted sample means and sample standard
deviations, not histogram bins or standard errors of the mean. This payload is independent of sampling and histogram retention. Joint output correlations are
not retained by the marginal surrogate. Access and transport reuse the stored
uncertainty-source identities.

Line-system construction checks exterior clearance on every realization. A draw
whose cables overlap throws `OverlapError` at its build.

Use one joint `Gridspace` builder for dependent internal dimensions. Derive
stack geometric boundaries and wire placements from its shared inputs. Repeating a Grid
as separate source arguments creates independent draws. The default normal
law has unbounded support: marginal means and standard uncertainties cannot
guarantee feasible geometry. `distribution=:uniform` instead has support
`nominal ± sqrt(3)*sigma` for each primitive descriptor. The builder must map
the complete joint support to feasible designs. Select this law only when the
uniform distribution is a justified statistical assumption. Changing the
marginal law neither repairs the geometry nor supplies a correlation.

`MonteCarlo(inner; options = (;), keywords...)` takes its execution controls in `options`
or as keywords of the same names. Both forms are normalized by
`computation_options(MonteCarlo, options)`, and a key supplied in both places is an error.

# Keywords

- `trials=nothing`: positive accepted-trial count, or DKW sizing when omitted.
- `confidence=0.95`: simultaneous empirical-CDF confidence [dimensionless], in `(0, 1)`.
- `cdf_tol=0.02`: maximum empirical-CDF deviation for DKW sizing [dimensionless], in `(0, 1)`.
- `distribution=:normal`: `:normal`, `:uniform`, a sampler function, or an extension-supported distribution.
- `seed=nothing`: nonnegative root seed representable as `UInt64`, or fresh randomness.
- `return_samples=false`: retain joint samples.
- `return_histograms=false`: retain marginal histogram densities.
- `bins=nothing`: positive histogram bin count, or automatic binning.
- `retain_details=false`: retain the computation details of each accepted trial.
- `max_failures=100`: positive maximum count of infeasible draws in a row for one trial.

The concrete type of `options` retains the sampler type.

$(TYPEDFIELDS)
"""
struct MonteCarlo{F <: AbstractFormulation, O <: ComputationOptions} <: AbstractTraversal
    "Formulation used for each sampled problem."
    inner::F
    "Normalized sampling and supplemental-output computation options."
    options::O

    MonteCarlo(inner::F, options::O) where {F <: AbstractFormulation, O <: ComputationOptions} =
        validate(new{F, O}(inner, options))
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
        ::Type{<:MonteCarlo},
        record::ComputationOptions
)::ComputationOptions
    options = record.data
    defaults = (
        trials = nothing, confidence = 0.95, cdf_tol = 0.02,
        distribution = :normal, seed = nothing,
        return_samples = false, return_histograms = false, bins = nothing,
        retain_details = false, max_failures = 100
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
        max_failures = Int(normalized.max_failures)
    )
end

function Base.NamedTuple(value::MonteCarlo)
    return (kind = :monte_carlo, inner = NamedTuple(value.inner), options = value.options.data)
end
