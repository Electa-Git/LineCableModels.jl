module Computation

export Combinatorial, LinearError, MonteCarlo, ParametricProblem
export ParametricResult, LinearErrorResult, MonteCarloResult
export CalculationManifest, ConfigurationFailure, SampleSummary, HistogramDensity, RLCG
export result, statistics, samples, histograms, uncertain_value, manifest

using Random
using SHA
using Statistics
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: PhaseDomain, basis, domain, frequencies
import ..DataModel
import ..Engine
import ..Grammar
import ..Grammar: compute
import ..ParametricBuilder
import ..PlotBuilder
import ..UnitHandler
using ..Grammar:
                 Combinatorial, LinearError, MonteCarlo, ParametricProblem,
                 ParametricResult, LinearErrorResult, MonteCarloResult,
                 CalculationManifest, ConfigurationFailure, SampleSummary,
                 HistogramDensity, RLCG,
                 Configuration, UncertainValue, AbsoluteUncertainty,
                 configurations, configuration_manifest, has_uncertainty,
                 materialize,
                 result, statistics, samples, histograms, uncertain_value, manifest

struct CableConstantsMaterializer end
function (::CableConstantsMaterializer)(design, separation, earth_resistivity)
    return Engine.CableConstantsProblem(design; separation, earth_resistivity)
end

function Engine.CableConstantsProblem(
        design::ParametricBuilder._AbstractDefinition{<:DataModel.CableDesign};
        separation = nothing,
        earth_resistivity = 100.0,
        combine::Symbol = :product
)
    return Grammar.Gridspace{Engine.CableConstantsProblem}(
        CableConstantsMaterializer(),
        (
            Grammar._gridspace_axis(design),
            Grammar._gridspace_axis(separation),
            Grammar._gridspace_axis(earth_resistivity)
        ),
        (:design, :separation, :earth_resistivity);
        combine
    )
end

function Grammar.SampleSummary(values::AbstractVector{<:Real})
    isempty(values) && throw(ArgumentError("cannot summarize an empty sample"))
    all(isfinite, values) || throw(ArgumentError("sample values must be finite"))
    sigma = length(values) == 1 ? zero(float(first(values))) : Statistics.std(values)
    promoted = promote(
        Statistics.mean(values),
        sigma,
        minimum(values),
        Statistics.quantile(values, 0.05),
        Statistics.quantile(values, 0.50),
        Statistics.quantile(values, 0.95),
        maximum(values)
    )
    min, q05, median, q95, max = promoted[3:7]
    min <= q05 <= median <= q95 <= max || throw(ArgumentError(
        "sample quantiles must be ordered between the minimum and maximum",
    ))
    promoted[2] >= zero(promoted[2]) || throw(ArgumentError(
        "sample standard deviation must be nonnegative",
    ))
    return SampleSummary(promoted..., length(values))
end

function Grammar.CalculationManifest(
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

mutable struct _ManifestState
    automatic_keys::IdDict{Any, Int}
end
_ManifestState() = _ManifestState(IdDict{Any, Int}())

function _manifest_tree(value, state::_ManifestState)
    if value isa Grammar.AutomaticGridKey
        label = get!(state.automatic_keys, value.token) do
            length(state.automatic_keys) + 1
        end
        return (automatic_grid = label,)
    elseif value isa Grammar.NamedGridKey
        return (named_grid = _manifest_tree(value.value, state),)
    elseif value isa Type
        return string(value)
    elseif value isa Function
        names = fieldnames(typeof(value))
        fields = map(name -> _manifest_tree(getfield(value, name), state), names)
        return (function_type = string(typeof(value)), fields = NamedTuple{names}(fields))
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
        return (
            type = string(typeof(value)),
            size = size(value),
            values = map(item -> _manifest_tree(item, state), Tuple(value))
        )
    elseif value isa Union{Nothing, Missing, Bool, Number, Symbol, AbstractString, Char}
        return value
    elseif isstructtype(typeof(value))
        names = fieldnames(typeof(value))
        fields = map(name -> _manifest_tree(getfield(value, name), state), names)
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
        ArgumentError, AssertionError, DimensionMismatch, DomainError
    }
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

function _details(; failures, manifest_value, samples_value = nothing,
        histograms_value = nothing, random_value = nothing)
    return Dict{Symbol, NamedTuple}(
        :failures => (values = failures,),
        :samples => (values = samples_value,),
        :histograms => (values = histograms_value,),
        :random => random_value === nothing ? (value = nothing,) : random_value,
        :manifest => (value = manifest_value,)
    )
end

function _manifest(problem, formulation, resolved, policy, options)
    return CalculationManifest(
        tuple(resolved...),
        _manifest_tree(problem.space),
        formulation,
        policy,
        options
    )
end

function _traverse(problem::ParametricProblem, formulation, result_type)
    values = nothing
    resolved = NamedTuple[]
    failures = ConfigurationFailure[]
    for (index, configuration) in enumerate(configurations(problem.space))
        parameterization = configuration_manifest(configuration)
        try
            materialized = materialize(configuration)
            values = _append_result(
                values,
                compute(materialized, formulation.inner; options = problem.options)
            )
            push!(resolved, parameterization)
        catch exception
            exception isa InterruptException && rethrow()
            (formulation.invalid === :error ||
             !_skippable_configuration_error(exception)) &&
                rethrow()
            push!(failures, _configuration_failure(index, configuration, exception))
        end
    end
    typed_values = values === nothing ? Grammar.AbstractProblemResult[] : values
    manifest_value = _manifest(
        problem, formulation.inner, resolved, formulation, problem.options)
    return result_type(
        formulation,
        typed_values,
        problem.space,
        _details(; failures, manifest_value)
    )
end

function compute(problem::ParametricProblem, formulation::Combinatorial)
    _traverse(problem, formulation, ParametricResult)
end

function compute(problem::ParametricProblem, formulation::LinearError)
    _traverse(problem, formulation, LinearErrorResult)
end

function _dkw_trials(observable_count::Integer, confidence::Real, cdf_tol::Real)
    observable_count > 0 || throw(ArgumentError("observable count must be positive"))
    alpha = 1 - confidence
    return ceil(Int, log(2 * observable_count / alpha) / (2 * cdf_tol^2))
end

_observable_count(::DataModel.CableConstants) = 3
function _observable_count(value::Engine.LineParameters)
    n = size(value.Z, 1)
    return 2 * n * (n + 1) * length(value.f)
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
    return HistogramDensity(edges, counts ./ (length(values) .* diff(edges)))
end

function _aggregate(draws::Vector{<:DataModel.CableConstants}, formulation::MonteCarlo)
    Rs = [value.R for value in draws]
    Ls = [value.L for value in draws]
    Cs = [value.C for value in draws]
    summaries = DataModel.CableConstants(SampleSummary(Rs), SampleSummary(Ls), SampleSummary(Cs))
    representation = DataModel.CableConstants(
        summaries.R.mean,
        summaries.L.mean,
        summaries.C.mean
    )
    retained = formulation.return_samples ? DataModel.CableConstants(Rs, Ls, Cs) : nothing
    hist = formulation.return_histograms ?
           DataModel.CableConstants(
        _histogram(Rs, formulation.bins),
        _histogram(Ls, formulation.bins),
        _histogram(Cs, formulation.bins)
    ) : nothing
    return (; representation, statistics = summaries, samples = retained, histograms = hist)
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

function _map_observables(function_value, sample_values::Array{<:Real, 4})
    output_size = size(sample_values)[1:3]
    indices = CartesianIndices(output_size)
    first_index = first(indices)
    first_value = function_value(collect(@view sample_values[first_index.I..., :]))
    output = Array{typeof(first_value)}(undef, output_size)
    output[first_index] = first_value
    for index in Iterators.drop(indices, 1)
        output[index] = function_value(collect(@view sample_values[index.I..., :]))
    end
    return output
end

function _aggregate(draws::Vector{<:Engine.LineParameters}, formulation::MonteCarlo)
    sample_values = _rlcg_samples(draws)
    summaries = RLCG((
        _map_observables(SampleSummary, values)
    for values in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
    )...)
    means = RLCG((
        map(summary -> summary.mean, values)
    for values in (summaries.R, summaries.L, summaries.C, summaries.G)
    )...)
    first_result = first(draws)
    angular = reshape(2π .* first_result.f, 1, 1, :)
    representation = Engine.LineParameters(
        domain(first_result),
        complex.(means.R, means.L .* angular),
        complex.(means.G, means.C .* angular),
        first_result.f;
        basis = basis(first_result)
    )
    hist = formulation.return_histograms ?
           RLCG((
        _map_observables(values -> _histogram(values, formulation.bins), samples)
    for samples in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
    )...) : nothing
    retained = formulation.return_samples ? sample_values : nothing
    return (; representation, statistics = summaries, samples = retained, histograms = hist)
end

function _monte_carlo(configuration, formulation::MonteCarlo, options, root_seed, seed)
    rng = Random.Xoshiro(seed)
    first_problem = rand(rng, configuration; distribution = formulation.distribution)
    first_result = compute(first_problem, formulation.inner; options)
    ntrials = something(
        formulation.trials,
        _dkw_trials(_observable_count(first_result), formulation.confidence, formulation.cdf_tol)
    )
    draws = Vector{typeof(first_result)}(undef, ntrials)
    draws[1] = first_result
    for trial in 2:ntrials
        realization = rand(rng, configuration; distribution = formulation.distribution)
        draws[trial] = compute(realization, formulation.inner; options)
    end
    return merge(_aggregate(draws, formulation), (; trials = ntrials, seed))
end

function compute(problem::ParametricProblem, formulation::MonteCarlo)
    values = nothing
    stats_values = nothing
    sample_values = Any[]
    histogram_values = Any[]
    resolved = NamedTuple[]
    failures = ConfigurationFailure[]
    seeds = UInt64[]
    trial_counts = Int[]
    root_seed = formulation.seed === nothing ? rand(Random.RandomDevice(), UInt64) :
                formulation.seed
    for (index, configuration) in enumerate(configurations(problem.space))
        parameterization = configuration_manifest(configuration)
        configuration_seed = root_seed ⊻ (UInt64(index - 1) * 0x9e3779b97f4a7c15)
        try
            aggregate = _monte_carlo(
                configuration,
                formulation,
                problem.options,
                root_seed,
                configuration_seed
            )
            values = _append_result(values, aggregate.representation)
            stats_values = _append_result(stats_values, aggregate.statistics)
            push!(sample_values, aggregate.samples)
            push!(histogram_values, aggregate.histograms)
            push!(seeds, aggregate.seed)
            push!(trial_counts, aggregate.trials)
            push!(resolved, parameterization)
        catch exception
            exception isa InterruptException && rethrow()
            (formulation.invalid === :error ||
             !_skippable_configuration_error(exception)) && rethrow()
            push!(failures, _configuration_failure(index, configuration, exception))
        end
    end
    typed_values = values === nothing ? Grammar.AbstractProblemResult[] : values
    typed_stats = stats_values === nothing ? Any[] : stats_values
    manifest_value = _manifest(
        problem, formulation.inner, resolved, formulation, problem.options)
    random_value = (
        root_seed = root_seed,
        configuration_seeds = seeds,
        trials = trial_counts,
        confidence = formulation.confidence,
        cdf_tol = formulation.cdf_tol,
        distribution = formulation.distribution
    )
    return MonteCarloResult(
        formulation,
        typed_values,
        problem.space,
        typed_stats,
        _details(
            ; failures,
            manifest_value,
            samples_value = formulation.return_samples ? sample_values : nothing,
            histograms_value = formulation.return_histograms ? histogram_values : nothing,
            random_value
        )
    )
end

include("dataframe.jl")
include("plotspecs.jl")

end
