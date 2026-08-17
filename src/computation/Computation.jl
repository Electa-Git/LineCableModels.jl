module Computation

export AbstractRunPolicy, FullParametric, MonteCarlo
export FullParametricResult, MonteCarloResult, CalculationManifest
export ConfigurationFailure, SampleSummary, HistogramPDF, RLCG
export result, statistics, samples, histograms, uncertain_value, manifest

using Random
using SHA
using Statistics
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: PhaseDomain, basis, domain, frequencies
import ..DataModel
import ..Engine
import ..Engine: compute!
import ..ParametricBuilder
import ..PlotBuilder
import ..UnitHandler
using ..ParametricBuilder:
                           AbstractSpec, Configuration, UncertainValue, AbsoluteUncertainty,
                           configurations, configuration_manifest, has_uncertainty,
                           materialize

struct CableConstantsMaterializer end
function (::CableConstantsMaterializer)(design, separation, earth_resistivity)
    return Engine.CableConstantsProblem(
        design;
        separation,
        earth_resistivity
    )
end

function Engine.CableConstantsProblem(
        design::AbstractSpec{<:DataModel.CableDesign};
        separation = nothing,
        earth_resistivity = 100.0,
        combine::Symbol = :product
)
    return ParametricBuilder.Gridspace{Engine.CableConstantsProblem}(
        CableConstantsMaterializer(),
        (
            ParametricBuilder._gridspace_axis(design),
            ParametricBuilder._gridspace_axis(separation),
            ParametricBuilder._gridspace_axis(earth_resistivity)
        ),
        (:design, :separation, :earth_resistivity);
        combine
    )
end

"""
$(TYPEDEF)

Supertype for policies that evaluate a formulation over parameter
configurations or stochastic realizations.
"""
abstract type AbstractRunPolicy end

"""
$(TYPEDEF)

Evaluate every configuration admitted by a recursively composed Gridspace.
Invalid configurations fail the traversal by default. `invalid=:skip` retains
an auditable [`ConfigurationFailure`](@ref) for each skipped configuration.

$(TYPEDFIELDS)
"""
struct FullParametric <: AbstractRunPolicy
    "Handling of invalid configurations: `:error` or `:skip`."
    invalid::Symbol

    function FullParametric(; invalid::Symbol = :error)
        invalid in (:error, :skip) || throw(ArgumentError(
            "invalid must be :error or :skip; got :$invalid",
        ))
        return new(invalid)
    end
end

"""
$(TYPEDEF)

Sample stochastic realizations conditional on each outer Gridspace
configuration. When `trials === nothing`, the trial count is selected with a
simultaneous Dvoretzky–Kiefer–Wolfowitz bound over the scalar observables.
`cdf_tol` is therefore an absolute empirical-CDF error bound, not a numerical
solver tolerance. With `M` scalar observables and confidence `1-α`, the policy
uses `ceil(log(2M/α) / (2*cdf_tol^2))`, so a union bound controls the maximum
empirical-CDF deviation across those marginals. It does not bound solver error,
mean error, or joint-distribution error.

$(TYPEDFIELDS)
"""
struct MonteCarlo{D, S} <: AbstractRunPolicy
    "Requested number of realizations, or `nothing` for DKW-based sizing."
    trials::Union{Nothing, Int}

    "Simultaneous marginal-CDF confidence level."
    confidence::Float64

    "Absolute empirical-CDF deviation used for automatic trial sizing."
    cdf_tol::Float64

    "Sampling distribution or compatible sampler."
    distribution::D

    "Root random seed, or `nothing` to obtain one from the system random device."
    seed::S

    "Whether complete observable samples are retained."
    return_samples::Bool

    "Whether marginal histogram densities are retained."
    return_histograms::Bool

    "Requested number of histogram bins, or `nothing` for automatic selection."
    bins::Union{Nothing, Int}

    "Handling of invalid configurations: `:error` or `:skip`."
    invalid::Symbol

    function MonteCarlo(;
            trials::Union{Nothing, Integer} = nothing,
            confidence::Real = 0.95,
            cdf_tol::Real = 0.02,
            distribution = :normal,
            seed::Union{Nothing, Integer} = nothing,
            return_samples::Bool = false,
            return_histograms::Bool = false,
            bins::Union{Nothing, Integer} = nothing,
            invalid::Symbol = :error
    )
        trials === nothing || trials > 0 ||
            throw(ArgumentError("trials must be positive"))
        0 < confidence < 1 ||
            throw(ArgumentError("confidence must lie between zero and one"))
        0 < cdf_tol < 1 ||
            throw(ArgumentError("cdf_tol must lie between zero and one"))
        bins === nothing || bins > 0 ||
            throw(ArgumentError("bins must be positive"))
        invalid in (:error, :skip) || throw(ArgumentError(
            "invalid must be :error or :skip; got :$invalid",
        ))
        distribution isa Symbol && distribution ∉ (:normal, :uniform) &&
            throw(ArgumentError(
                "unsupported distribution :$distribution; expected :normal, :uniform, a sampler function, or an extension-supported distribution",
            ))
        actual_seed = seed === nothing ? nothing : UInt64(seed)
        return new{typeof(distribution), typeof(actual_seed)}(
            trials === nothing ? nothing : Int(trials),
            Float64(confidence),
            Float64(cdf_tol),
            distribution,
            actual_seed,
            return_samples,
            return_histograms,
            bins === nothing ? nothing : Int(bins),
            invalid
        )
    end
end

"""
$(TYPEDEF)

Store descriptive statistics for one scalar Monte Carlo observable.

$(TYPEDFIELDS)
"""
struct SampleSummary{T <: Real}
    "Sample mean, in the physical unit of the observable."
    mean::T

    "Sample standard deviation, in the physical unit of the observable."
    std::T

    "Minimum sampled value, in the physical unit of the observable."
    min::T

    "Fifth sample percentile, in the physical unit of the observable."
    q05::T

    "Median sampled value, in the physical unit of the observable."
    q50::T

    "Ninety-fifth sample percentile, in the physical unit of the observable."
    q95::T

    "Maximum sampled value, in the physical unit of the observable."
    max::T

    function SampleSummary(
            mean::T,
            std::T,
            min::T,
            q05::T,
            q50::T,
            q95::T,
            max::T
    ) where {T <: Real}
        all(isfinite, (mean, std, min, q05, q50, q95, max)) ||
            throw(ArgumentError("sample summary values must be finite"))
        std >= zero(std) || throw(ArgumentError(
            "sample standard deviation must be nonnegative",
        ))
        min <= q05 <= q50 <= q95 <= max || throw(ArgumentError(
            "sample quantiles must be ordered between the minimum and maximum",
        ))
        min <= mean <= max || throw(ArgumentError(
            "sample mean must lie between the minimum and maximum",
        ))
        return new{T}(mean, std, min, q05, q50, q95, max)
    end
end

"""
$(TYPEDSIGNATURES)

Calculate a [`SampleSummary`](@ref) from a nonempty vector of finite values.
"""
function SampleSummary(values::AbstractVector{<:Real})
    isempty(values) && throw(ArgumentError("cannot summarize an empty sample"))
    all(isfinite, values) || throw(ArgumentError("sample values must be finite"))
    sigma = length(values) == 1 ? zero(float(first(values))) : Statistics.std(values)
    promoted = promote(
        Statistics.mean(values), sigma, minimum(values),
        Statistics.quantile(values, 0.05), Statistics.quantile(values, 0.50),
        Statistics.quantile(values, 0.95), maximum(values)
    )
    return SampleSummary(promoted...)
end

"""
$(TYPEDEF)

Represent a normalized piecewise-constant marginal probability density.

$(TYPEDFIELDS)
"""
struct HistogramPDF{T <: AbstractFloat}
    "Bin edges in the physical unit of the observable."
    edges::Vector{T}

    "Probability densities, in the reciprocal unit of the observable."
    density::Vector{T}

    function HistogramPDF(edges::Vector{T}, density::Vector{T}) where {T <: AbstractFloat}
        length(edges) == length(density) + 1 || throw(ArgumentError(
            "histogram edges must contain one more value than density",
        ))
        isempty(density) && throw(ArgumentError("histogram density cannot be empty"))
        all(isfinite, edges) ||
            throw(ArgumentError("histogram edges must be finite"))
        widths = diff(edges)
        all(>(zero(T)), widths) ||
            throw(ArgumentError("histogram edges must be strictly increasing"))
        all(x -> isfinite(x) && x >= zero(T), density) ||
            throw(ArgumentError("histogram density must be finite and nonnegative"))
        area = sum(density .* widths)
        area > zero(T) || throw(ArgumentError("histogram density must have positive area"))
        return new{T}(copy(edges), density ./ area)
    end
end

"""
$(TYPEDSIGNATURES)

Construct and normalize a piecewise-constant histogram density.
"""
function HistogramPDF(
        edges::AbstractVector{<:Real},
        density::AbstractVector{<:Real}
)
    T = promote_type(float(eltype(edges)), float(eltype(density)))
    isconcretetype(T) || (T = Float64)
    return HistogramPDF(Vector{T}(edges), Vector{T}(density))
end

"""
$(TYPEDEF)

Group resistance, inductance, capacitance, and conductance representations.

$(TYPEDFIELDS)
"""
struct RLCG{T}
    "Resistance representation `\\[Ω/m\\]` or `\\[Ω\\]`, according to result basis."
    R::T

    "Inductance representation `\\[H/m\\]` or `\\[H\\]`, according to result basis."
    L::T

    "Capacitance representation `\\[F/m\\]` or `\\[F\\]`, according to result basis."
    C::T

    "Conductance representation `\\[S/m\\]` or `\\[S\\]`, according to result basis."
    G::T
end

"""
$(TYPEDEF)

Record the inputs that define one calculation or complete Gridspace traversal
and their stable digest.

$(TYPEDFIELDS)
"""
struct CalculationManifest{R, P, F, S, E, O}
    "Resolved Gridspace parameterization."
    resolved_parameterization::R

    "Materialized problem assumptions before evaluation."
    problem_assumptions::P

    "Formulation and formulation options."
    formulation::F

    "Numerical solver identifier."
    solver::S

    "Configuration or realization evaluation policy."
    execution_policy::E

    "Computation options supplied to `compute!`."
    calculation_options::O

    "SHA-256 digest of the stable manifest serialization."
    hash::String
end

"""
$(TYPEDSIGNATURES)

Create a reproducible calculation manifest and its SHA-256 digest.
"""
function CalculationManifest(
        resolved_parameterization,
        problem_assumptions,
        formulation,
        execution_policy,
        calculation_options
)
    solver = string(typeof(formulation))
    payload = (
        resolved_parameterization = resolved_parameterization,
        problem_assumptions = problem_assumptions,
        formulation = _manifest_tree(formulation),
        solver,
        execution_policy = _manifest_tree(execution_policy),
        calculation_options = _manifest_tree(calculation_options)
    )
    digest = bytes2hex(sha256(_stable_bytes(payload)))
    return CalculationManifest(
        payload.resolved_parameterization,
        payload.problem_assumptions,
        payload.formulation,
        payload.solver,
        payload.execution_policy,
        payload.calculation_options,
        digest
    )
end

"""
$(TYPEDEF)

Record an invalid configuration skipped by explicit request.

$(TYPEDFIELDS)
"""
struct ConfigurationFailure{C}
    "One-based position of the configuration in the traversal."
    index::Int

    "Resolved parameter values for the configuration."
    configuration::C

    "Qualified exception type name."
    exception_type::String

    "Rendered exception message."
    message::String
end

"""
$(TYPEDEF)

Store the results and parameter values from a complete Gridspace traversal.

$(TYPEDFIELDS)
"""
struct FullParametricResult{T, C, F, M} <: AbstractVector{T}
    "Primitive result for each successfully evaluated configuration."
    values::Vector{T}

    "Resolved parameterization corresponding to each value."
    configurations::Vector{C}

    "Invalid configurations skipped by explicit request."
    failures::Vector{F}

    "Manifest describing the complete traversal."
    manifest::M
end

Base.size(value::FullParametricResult) = size(value.values)
Base.length(value::FullParametricResult) = length(value.values)
Base.getindex(value::FullParametricResult, index::Integer) = value.values[index]
Base.iterate(value::FullParametricResult, state...) = iterate(value.values, state...)
Base.IndexStyle(::Type{<:FullParametricResult}) = IndexLinear()

"""
$(TYPEDEF)

Store one Monte Carlo analysis of primitive result type `T`.

$(TYPEDFIELDS)
"""
struct MonteCarloResult{T, S, Sa, H, U, D, M}
    "Complete primitive result reconstructed from the sample means."
    representation::T

    "Marginal descriptive statistics in the shape of `T`."
    statistics::S

    "Retained realization values, or `nothing`."
    samples::Sa

    "Retained marginal histogram densities, or `nothing`."
    histograms::H

    "Standard-uncertainty representation in the shape of `T`."
    uncertain::U

    "Number of stochastic realizations."
    trials::Int

    "Simultaneous marginal-CDF confidence level used for automatic sizing."
    confidence::Float64

    "Absolute empirical-CDF deviation used for automatic sizing."
    cdf_tol::Float64

    "Distribution or sampler used to draw realizations."
    distribution::D

    "Random seed used for this configuration."
    seed::UInt64

    "Manifest describing this Monte Carlo analysis."
    manifest::M
end

"""
$(TYPEDSIGNATURES)

Return the complete primitive representation stored in a computation result.
"""
result(value::FullParametricResult) = value.values
result(value::MonteCarloResult) = value.representation

"""
$(TYPEDSIGNATURES)

Return the marginal descriptive statistics of a Monte Carlo result.
"""
statistics(value::MonteCarloResult) = value.statistics

"""
$(TYPEDSIGNATURES)

Return retained realizations, or `nothing` when they were not requested.
"""
samples(value::MonteCarloResult) = value.samples

"""
$(TYPEDSIGNATURES)

Return retained marginal histogram densities, or `nothing` when they were not
requested.
"""
histograms(value::MonteCarloResult) = value.histograms

"""
$(TYPEDSIGNATURES)

Return the standard-uncertainty representation of a Monte Carlo result.
"""
uncertain_value(value::MonteCarloResult) = value.uncertain

"""
$(TYPEDSIGNATURES)

Return the calculation manifest associated with a result.
"""
manifest(value::Union{FullParametricResult, MonteCarloResult}) = value.manifest

# Stable serialization deliberately avoids Julia's session-randomized hash
# and object addresses. Automatic Grid coupling identities are assigned stable
# first-occurrence labels during structural traversal.
mutable struct _ManifestState
    automatic_keys::IdDict{Any, Int}
end
_ManifestState() = _ManifestState(IdDict{Any, Int}())

function _manifest_tree(value, state::_ManifestState)
    if value isa ParametricBuilder.AutomaticGridKey
        label = get!(state.automatic_keys, value.token) do
            length(state.automatic_keys) + 1
        end
        return (automatic_grid = label,)
    elseif value isa ParametricBuilder.NamedGridKey
        return (named_grid = _manifest_tree(value.value, state),)
    elseif value isa Type
        return string(value)
    elseif value isa Function
        names = fieldnames(typeof(value))
        fields = map(
            name -> _manifest_tree(getfield(value, name), state),
            names
        )
        return (
            function_type = string(typeof(value)),
            fields = NamedTuple{names}(fields)
        )
    elseif value isa NamedTuple
        return NamedTuple{keys(value)}(map(item -> _manifest_tree(item, state), values(value)))
    elseif value isa AbstractDict
        entries = [(_manifest_tree(key, state), _manifest_tree(item, state))
                   for (key, item) in pairs(value)]
        sort!(entries; by = entry -> String(_stable_bytes(first(entry))))
        return (type = string(typeof(value)), entries = tuple(entries...))
    elseif value isa AbstractSet
        entries = [_manifest_tree(item, state) for item in value]
        sort!(entries; by = entry -> String(_stable_bytes(entry)))
        return (type = string(typeof(value)), entries = tuple(entries...))
    elseif value isa Tuple
        return map(item -> _manifest_tree(item, state), value)
    elseif value isa AbstractArray
        return (type = string(typeof(value)), size = size(value),
            values = map(
                item -> _manifest_tree(item, state), Tuple(value)
            ))
    elseif value isa Union{Nothing, Missing, Bool, Number, Symbol, AbstractString, Char}
        return value
    elseif isstructtype(typeof(value))
        names = fieldnames(typeof(value))
        fields = map(
            name -> _manifest_tree(getfield(value, name), state),
            names
        )
        return (type = string(typeof(value)), fields = NamedTuple{names}(fields))
    end
    return (type = string(typeof(value)), representation = repr(value))
end

_manifest_tree(value) = _manifest_tree(value, _ManifestState())

function _stable_bytes(value)
    io = IOBuffer()
    _write_stable(io, value)
    return take!(io)
end

function _write_stable(io::IO, value)
    if value isa NamedTuple
        print(io, "named{")
        for key in keys(value)
            print(io, repr(key), '=')
            _write_stable(io, getproperty(value, key))
            print(io, ';')
        end
        print(io, '}')
    elseif value isa Tuple
        print(io, "tuple[")
        for item in value
            _write_stable(io, item)
            print(io, ';')
        end
        print(io, ']')
    elseif value isa AbstractArray
        print(io, "array", repr(size(value)), '[')
        for item in value
            _write_stable(io, item)
            print(io, ';')
        end
        print(io, ']')
    elseif value isa AbstractFloat
        print(io, string(typeof(value)), ':', bitstring(value))
    elseif value isa Complex
        print(io, "complex(")
        _write_stable(io, real(value))
        print(io, ',')
        _write_stable(io, imag(value))
        print(io, ')')
    elseif value isa Number
        print(io, string(typeof(value)), ':', repr(value))
    else
        print(io, string(typeof(value)), ':', repr(value))
    end
    return io
end

function _configuration_failure(index, configuration, exception)
    return ConfigurationFailure(
        index,
        configuration_manifest(configuration),
        string(typeof(exception)),
        sprint(showerror, exception)
    )
end

function _skippable_configuration_error(exception)
    exception isa Union{
        ArgumentError,
        AssertionError,
        DimensionMismatch,
        DomainError
    }
end

function _full_result(values, resolved, failures, problem, formulation, run, options)
    typed_values = values === nothing ? Union{}[] : values
    calculation_manifest = CalculationManifest(
        tuple(resolved...),
        _manifest_tree(problem),
        formulation,
        run,
        options
    )
    return FullParametricResult(
        typed_values,
        resolved,
        failures,
        calculation_manifest
    )
end

_append_result(::Nothing, value) = typeof(value)[value]
function _append_result(values::Vector{T}, value) where {T}
    value isa T && (push!(values, value); return values)
    W = typejoin(T, typeof(value))
    W === Any && throw(ArgumentError(
        "configuration results do not share one primitive result grammar",
    ))
    widened = Vector{W}(undef, length(values) + 1)
    copyto!(widened, values)
    widened[end] = value
    return widened
end

"""
$(TYPEDSIGNATURES)

Evaluate a formulation over every configuration admitted by a problem
specification.

Deterministic specifications use [`FullParametric`](@ref) when `run` is
omitted. Specifications containing uncertainty require an explicit
[`FullParametric`](@ref) or [`MonteCarlo`](@ref) policy.
"""
function compute!(
        problem::AbstractSpec{<:Engine.ProblemDefinition},
        formulation::Engine.AbstractFormulation;
        run = nothing,
        options = Engine.ComputeOptions()
)
    selected_run = if run === nothing
        has_uncertainty(problem) && throw(ArgumentError(
            "uncertain Gridspaces require explicit run=FullParametric() or run=MonteCarlo(...) execution intent",
        ))
        FullParametric()
    else
        run
    end
    selected_run isa AbstractRunPolicy || throw(ArgumentError(
        "run must be FullParametric() or MonteCarlo(...)",
    ))
    normalized_options = Engine.compute_options(options)
    return _compute_gridspace(problem, formulation, selected_run, normalized_options)
end

function _compute_gridspace(problem, formulation, run::FullParametric, options)
    values = nothing
    resolved = NamedTuple[]
    failures = ConfigurationFailure[]
    for (index, configuration) in enumerate(configurations(problem))
        parameterization = configuration_manifest(configuration)
        try
            materialized = materialize(configuration)
            values = _append_result(
                values,
                Engine.compute!(materialized, formulation; options)
            )
            push!(resolved, parameterization)
        catch exception
            exception isa InterruptException && rethrow()
            (run.invalid === :error || !_skippable_configuration_error(exception)) &&
                rethrow()
            push!(failures, _configuration_failure(index, configuration, exception))
        end
    end
    return _full_result(values, resolved, failures, problem, formulation, run, options)
end

function _compute_gridspace(problem, formulation, run::MonteCarlo, options)
    values = nothing
    resolved = NamedTuple[]
    failures = ConfigurationFailure[]
    root_seed = run.seed === nothing ? rand(Random.RandomDevice(), UInt64) : run.seed
    for (index, configuration) in enumerate(configurations(problem))
        parameterization = configuration_manifest(configuration)
        configuration_seed = root_seed ⊻ (UInt64(index - 1) * 0x9e3779b97f4a7c15)
        try
            values = _append_result(
                values,
                _monte_carlo(
                    problem,
                    configuration,
                    parameterization,
                    formulation,
                    run,
                    options,
                    root_seed,
                    configuration_seed
                )
            )
            push!(resolved, parameterization)
        catch exception
            exception isa InterruptException && rethrow()
            (run.invalid === :error || !_skippable_configuration_error(exception)) &&
                rethrow()
            push!(failures, _configuration_failure(index, configuration, exception))
        end
    end
    manifest_run = (
        policy = _manifest_tree(run),
        actual_root_seed = root_seed
    )
    values === nothing && return _full_result(
        values, resolved, failures, problem, formulation, manifest_run, options
    )
    if length(values) == 1 && isempty(failures)
        return only(values)
    end
    return _full_result(
        values,
        resolved,
        failures,
        problem,
        formulation,
        manifest_run,
        options
    )
end

function _dkw_trials(observables::Integer, confidence::Real, cdf_tol::Real)
    observables > 0 || throw(ArgumentError("observable count must be positive"))
    alpha = 1 - confidence
    return ceil(Int, log(2 * observables / alpha) / (2 * cdf_tol^2))
end

_observable_count(::DataModel.CableConstants) = 3
function _observable_count(value::Engine.LineParameters)
    n = size(value.Z, 1)
    return 2 * n * (n + 1) * length(value.f)
end

function _monte_carlo(
        problem,
        configuration,
        parameterization,
        formulation,
        run,
        options,
        root_seed,
        seed
)
    rng = Random.Xoshiro(seed)
    first_problem = rand(rng, configuration; distribution = run.distribution)
    first_result = Engine.compute!(first_problem, formulation; options)
    ntrials = something(
        run.trials,
        _dkw_trials(_observable_count(first_result), run.confidence, run.cdf_tol)
    )
    draws = Vector{typeof(first_result)}(undef, ntrials)
    draws[1] = first_result
    for trial in 2:ntrials
        realization = rand(rng, configuration; distribution = run.distribution)
        draws[trial] = Engine.compute!(realization, formulation; options)
    end
    aggregated = _aggregate(draws, run)
    calculation_manifest = CalculationManifest(
        parameterization,
        _manifest_tree(problem),
        formulation,
        merge(
            _manifest_tree(run),
            (
                actual_trials = ntrials,
                actual_root_seed = root_seed,
                actual_seed = seed
            )
        ),
        options
    )
    return MonteCarloResult(
        aggregated.representation,
        aggregated.statistics,
        aggregated.samples,
        aggregated.histograms,
        UncertainValue(
            aggregated.representation,
            aggregated.uncertainty,
            AbsoluteUncertainty()
        ),
        ntrials,
        run.confidence,
        run.cdf_tol,
        run.distribution,
        seed,
        calculation_manifest
    )
end

function _histogram(values::AbstractVector{<:Real}, bins::Union{Nothing, Int})
    lo, hi = extrema(values)
    if lo == hi
        width = max(abs(float(lo)) * sqrt(eps(Float64)), sqrt(eps(Float64)))
        edges = [float(lo) - width, float(hi) + width]
    else
        count = something(bins, max(1, ceil(Int, sqrt(length(values)))))
        edges = collect(range(float(lo), float(hi); length = count + 1))
    end
    counts = zeros(Float64, length(edges) - 1)
    for value in values
        index = value == last(edges) ? length(counts) :
                clamp(searchsortedlast(edges, value), 1, length(counts))
        counts[index] += 1
    end
    density = counts ./ (length(values) .* diff(edges))
    return HistogramPDF(edges, density)
end

function _aggregate(draws::Vector{<:DataModel.CableConstants}, run)
    Rs = [value.R for value in draws]
    Ls = [value.L for value in draws]
    Cs = [value.C for value in draws]
    summaries = DataModel.CableConstants(
        SampleSummary(Rs), SampleSummary(Ls), SampleSummary(Cs)
    )
    representation = DataModel.CableConstants(
        summaries.R.mean, summaries.L.mean, summaries.C.mean
    )
    uncertainty = DataModel.CableConstants(
        summaries.R.std, summaries.L.std, summaries.C.std
    )
    retained = run.return_samples ? DataModel.CableConstants(Rs, Ls, Cs) : nothing
    hist = run.return_histograms ?
           DataModel.CableConstants(
        _histogram(Rs, run.bins),
        _histogram(Ls, run.bins),
        _histogram(Cs, run.bins)
    ) : nothing
    return (representation, statistics = summaries, samples = retained,
        histograms = hist, uncertainty)
end

function _rlcg_samples(draws::Vector{<:Engine.LineParameters})
    first_result = first(draws)
    dimensions = (size(first_result.Z)..., length(draws))
    Rs = Array{Float64}(undef, dimensions)
    Ls = similar(Rs)
    Cs = similar(Rs)
    Gs = similar(Rs)
    angular = reshape(2π .* first_result.f, 1, 1, :)
    for (trial, value) in enumerate(draws)
        size(value.Z) == size(first_result.Z) || throw(DimensionMismatch(
            "Monte Carlo realizations produced incompatible impedance dimensions",
        ))
        value.f == first_result.f || throw(DimensionMismatch(
            "Monte Carlo realizations produced incompatible frequency axes",
        ))
        @views Rs[:, :, :, trial] .= real.(value.Z.values)
        @views Ls[:, :, :, trial] .= imag.(value.Z.values) ./ angular
        @views Gs[:, :, :, trial] .= real.(value.Y.values)
        @views Cs[:, :, :, trial] .= imag.(value.Y.values) ./ angular
    end
    return RLCG(Rs, Ls, Cs, Gs)
end

function _map_observables(function_value, samples::Array{<:Real, 4})
    output_size = size(samples)[1:3]
    indices = CartesianIndices(output_size)
    first_index = first(indices)
    first_value = function_value(collect(@view samples[first_index.I..., :]))
    output = Array{typeof(first_value)}(undef, output_size)
    output[first_index] = first_value
    for index in Iterators.drop(indices, 1)
        output[index] = function_value(collect(@view samples[index.I..., :]))
    end
    return output
end

function _aggregate(draws::Vector{<:Engine.LineParameters}, run)
    sample_values = _rlcg_samples(draws)
    summaries = RLCG((
        _map_observables(SampleSummary, values)
    for values in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
    )...)

    means = RLCG((map(summary -> summary.mean, values)
    for values in (summaries.R, summaries.L, summaries.C, summaries.G))...)
    deviations = RLCG((map(summary -> summary.std, values)
    for values in (summaries.R, summaries.L, summaries.C, summaries.G))...)
    first_result = first(draws)
    angular = reshape(2π .* first_result.f, 1, 1, :)
    mean_Z = complex.(means.R, means.L .* angular)
    mean_Y = complex.(means.G, means.C .* angular)
    representation = Engine.LineParameters(
        domain(first_result), mean_Z, mean_Y, first_result.f;
        basis = basis(first_result)
    )
    hist = run.return_histograms ?
           RLCG((
        _map_observables(
            values -> _histogram(values, run.bins),
            samples
        )
    for samples in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
    )...) : nothing
    retained = run.return_samples ? sample_values : nothing
    return (representation, statistics = summaries, samples = retained,
        histograms = hist, uncertainty = deviations)
end

include("dataframe.jl")
include("plotspecs.jl")

end
