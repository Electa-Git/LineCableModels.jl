"""
$(TYPEDEF)

Select direct linear uncertainty propagation with `inner`.

Shared inputs must be supplied once to a joint `Gridspace` builder, which
derives the dependent geometry. Propagation retains those correlations and
differentiates continuous geometry at the nominal design, including bounded
power-cell compaction. Strand counts and clipping topology are nominal discrete
choices. A smooth linear model cannot describe a study across such a transition.
First-order moments need not equal nonlinear Monte Carlo moments.

The study has no options of its own: the computation options of its run are those of the
parametric call, `compute(problem, study; options)`. A point fails as a point of a
[`Combinatorial`](@ref) study does, and the result records it in its `failures`.

$(TYPEDFIELDS)
"""
struct LinearError{F <: AbstractFormulation} <: AbstractTraversal
    "Formulation used for each materialized problem."
    inner::F

    LinearError(inner::F) where {F <: AbstractFormulation} = validate(new{F}(inner))
end

# The record of linear propagation contains its inner formulation and an empty set of options.
function Base.NamedTuple(value::LinearError)
    return (kind = :linear_error, inner = NamedTuple(value.inner), options = (;))
end

"""Identify first-order uncertainty propagation without inspecting its inner owner."""
description(::Type{<:LinearError}; compact::Bool=false) = "LEP"
description(::LinearError; compact::Bool=false) = description(LinearError;compact)
formula_id(::Type{<:LinearError}) = :LinearError
formula_id(::LinearError) = :LinearError
function Base.pairs(value::LinearError; quantity = nothing)
    pairs(LinearError, (inner = value.inner, options = (;)); quantity)
end

"""
$(TYPEDEF)

Select conditional Monte Carlo propagation over a
[`ParametricProblem`](@ref). Randomness is local and reproducible when `seed`
is supplied.

A study computes its points through the runner, `traverse`. Before the first draw, it
builds the nominal problem of each point from the means of its uncertain arguments, and the
constructor of the problem validates it. Every point then computes its trials from its
uncertain arguments with the backend of the inner formulation and returns its own
`MonteCarloResult`, and the result of the study joins these results. The run logs the
timing of each point, and its progress between the points and between the trials of a
point. After each computed point, `on_result(problem, index, point)` receives the nominal
problem, the position of the point in the space and the `MonteCarloResult` of that point.
`verbosity`, `on_result` and `retain_details`, which retains the computation details of each
accepted trial, are the computation options of the run, passed in the parametric call
`compute(problem, study; options)`.

A draw is infeasible when its build throws
[`OverlapError`](@ref LineCableModels.DataModel.OverlapError), or when its build or
computation throws `DomainError`: [`on_error`](@ref LineCableModels.Commons.on_error)
returns these exceptions for a `MonteCarlo` study. The study discards an infeasible draw as
if it was never drawn, records nothing about it and draws again. The output distribution is
conditional on feasible draws.

After `max_failures` infeasible draws in a row for one trial, the point fails with a
`DomainError`: the distribution of the point has no feasible draw. Its value is the
uncertain arguments of the point, and its message gives that count and the last exception.
When `on_error` returns the exception of a point for the study, as it does for this
`DomainError` and for an `OverlapError` of a nominal problem, the result records the point in
its `failures`, and the study continues. Any other exception stops the study.

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

`MonteCarlo(inner; options = (;), keywords...)` takes the options of the study in `options`
or as keywords of the same names, and a key supplied in both places is an error. The
constructor `MonteCarlo(inner, options)` normalizes them with
`formulation_options(MonteCarlo, options)`, and `formulation_options(study)` returns them.
They are scalars: a `Grid` among them makes a parameter space of option records, which no
method of `MonteCarlo` accepts.

# Keywords

- `trials=nothing`: positive accepted-trial count, or DKW sizing when omitted.
- `confidence=0.95`: simultaneous empirical-CDF confidence [dimensionless], in `(0, 1)`.
- `cdf_tol=0.02`: maximum empirical-CDF deviation for DKW sizing [dimensionless], in `(0, 1)`.
- `distribution=:normal`: `:normal`, `:uniform`, a sampler function, or an extension-supported distribution.
- `seed=nothing`: nonnegative root seed representable as `UInt64`, or fresh randomness.
- `return_samples=false`: retain joint samples.
- `return_histograms=false`: retain marginal histogram densities.
- `bins=nothing`: positive histogram bin count, or automatic binning.
- `max_failures=100`: positive maximum count of infeasible draws in a row for one trial.

The concrete type of `options` retains the sampler type.

$(TYPEDFIELDS)
"""
struct MonteCarlo{F <: AbstractFormulation, O <: FormulationOptions} <: AbstractTraversal
    "Formulation used for each sampled problem."
    inner::F
    "The normalized options of the study: its sampling and what its result retains."
    options::O

    function MonteCarlo(inner::F, record::FormulationOptions) where {F <: AbstractFormulation}
        options = formulation_options(MonteCarlo, record)
        return validate(new{F, typeof(options)}(inner, options))
    end
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
function Base.pairs(value::MonteCarlo; quantity = nothing)
    pairs(MonteCarlo, (inner = value.inner, options = value.options.data); quantity)
end
function Base.pairs(owner::Type{<:Union{MonteCarlo, LinearError}}, retained::NamedTuple; quantity = nothing)
    entries=Pair{Tuple, Any}[(owner, ()) => (owner => retained.options)]
    append!(entries, pairs(retained.inner; quantity))
    return entries
end
description(::Type{<:Union{MonteCarlo, LinearError}}, ::Val{:inner}; compact::Bool=false) = "inner method"

function MonteCarlo(inner; options::NamedTuple = (;), keywords...)
    duplicates = filter(key -> haskey(options, key), keys(keywords))
    isempty(duplicates) || throw(ArgumentError(
        "MonteCarlo options supplied both as keywords and in options: $(collect(duplicates))",
    ))
    return MonteCarlo(inner, FormulationOptions(merge(options, values(keywords))))
end

# Normalize the options of a Monte Carlo study. The record stores the values themselves, so a
# sampler function stays the same object.
function formulation_options(
        ::Type{<:MonteCarlo},
        record::FormulationOptions
)::FormulationOptions
    options = record.data
    defaults = (
        trials = nothing, confidence = 0.95, cdf_tol = 0.02,
        distribution = :normal, seed = nothing,
        return_samples = false, return_histograms = false, bins = nothing,
        max_failures = 100
    )
    unknown = filter(key -> key ∉ keys(defaults), keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown MonteCarlo options: $(sort!(collect(unknown))). The computation options " *
        "of the run are those of compute(problem, study; options)",
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
    for key in (:return_samples, :return_histograms)
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
    data = (
        trials = normalized.trials === nothing ? nothing : Int(normalized.trials),
        confidence = Float64(normalized.confidence),
        cdf_tol = Float64(normalized.cdf_tol),
        distribution = distribution,
        seed = seed === nothing ? nothing : UInt64(seed),
        return_samples = normalized.return_samples,
        return_histograms = normalized.return_histograms,
        bins = normalized.bins === nothing ? nothing : Int(normalized.bins),
        max_failures = Int(normalized.max_failures)
    )
    return FormulationOptions{typeof(data)}(data)
end
formulation_options(study::MonteCarlo) = study.options

function Base.NamedTuple(value::MonteCarlo)
    return (kind = :monte_carlo, inner = NamedTuple(value.inner), options = value.options.data)
end
