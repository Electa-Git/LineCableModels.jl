function _dkw_trials(observable_count::Integer, confidence::Real, cdf_tol::Real)
    observable_count > 0 || throw(ArgumentError("observable count must be positive"))
    alpha = 1 - confidence
    return ceil(Int, log(2 * observable_count / alpha) / (2 * cdf_tol^2))
end

_observable_count(value::Engine.CableConstants) = 4length(value)

function _observable_count(value::Engine.LineParameters)
    n = size(observe(value, Engine.Z), 1)
    return 2 * n * (n + 1) * length(observe(value, Engine.frequencies))
end

"""
$(TYPEDSIGNATURES)

Extend `buffers` with the `R`, `L`, `C` and `G` samples of `plan.trials` realizations,
one row per assembly of `input`, the first accepted result. The samples keep the scalar
type of that result.
"""
function initialize_buffers(::MonteCarlo, ::Type, input::Engine.CableConstants, plan,
        buffers)
    T = eltype(observe(input, R))
    return merge(buffers, (
        R = Matrix{T}(undef, length(input), plan.trials),
        L = Matrix{T}(undef, length(input), plan.trials),
        C = Matrix{T}(undef, length(input), plan.trials),
        G = Matrix{T}(undef, length(input), plan.trials)
    ))
end

function _record_sample!(
        buffers::NamedTuple{(:R, :L, :C, :G)},
        value::Engine.CableConstants,
        trial::Int,
        expected_axis
)
    value.terminals == expected_axis.terminals || throw(DimensionMismatch(
        "Monte Carlo cable-constant realizations produced incompatible assembly labels",
    ))
    value.frequency == expected_axis.frequency || throw(DimensionMismatch(
        "Monte Carlo cable-constant realizations produced incompatible frequencies",
    ))
    buffers.R[:, trial] .= observe(value, R)
    buffers.L[:, trial] .= observe(value, L)
    buffers.C[:, trial] .= observe(value, C)
    buffers.G[:, trial] .= observe(value, Engine.G)
    return buffers
end

function _sample_axis(value::Engine.CableConstants)
    (
        terminals = copy(value.terminals),
        frequency = value.frequency
    )
end
_sample_axis(value::Engine.LineParameters) = observe(value, Engine.frequencies)

function _cable_summaries(values::AbstractMatrix)
    return [SampleSummary(collect(@view values[assembly, :]))
            for assembly in axes(values, 1)]
end

function _sample_shunt_model(value::Union{Engine.LineParameters,Engine.CableConstants})
    model=get(details(value).data,:shunt_model,nothing)
    model===nothing && return ComputationDetails()
    return ComputationDetails(NamedTuple{(:shunt_model,),Tuple{NamedTuple}}(((
        requested=model.requested,
        domains=[(design=d.design,terminals=collect(d.terminals),reason=d.reason)
            for d in model.domains]),)))
end

function _aggregate(
        sample_values::NamedTuple{(:R, :L, :C, :G)},
        first_result::Engine.CableConstants,
        formulation::MonteCarlo
)
    summaries = (
        R = _cable_summaries(sample_values.R),
        L = _cable_summaries(sample_values.L),
        C = _cable_summaries(sample_values.C),
        G = _cable_summaries(sample_values.G)
    )
    representation = materialize(first_result, summaries; details=_sample_shunt_model(first_result))
    retained = formulation.options.data.return_samples ? sample_values : nothing
    hist = formulation.options.data.return_histograms ?
           (
        R = [HistogramDensity(collect(@view sample_values.R[assembly, :]);
                 bins = formulation.options.data.bins) for assembly in axes(sample_values.R, 1)],
        L = [HistogramDensity(collect(@view sample_values.L[assembly, :]);
                 bins = formulation.options.data.bins) for assembly in axes(sample_values.L, 1)],
        C = [HistogramDensity(collect(@view sample_values.C[assembly, :]);
                 bins = formulation.options.data.bins) for assembly in axes(sample_values.C, 1)],
        G = [HistogramDensity(collect(@view sample_values.G[assembly, :]);
                 bins = formulation.options.data.bins) for assembly in axes(sample_values.G, 1)]
    ) : nothing
    return (; representation, statistics = summaries, samples = retained, histograms = hist)
end

"""
$(TYPEDSIGNATURES)

Extend `buffers` with the `R`, `L`, `C` and `G` samples of `plan.trials` realizations, in
the phase and frequency layout of `input`, the first accepted result, with scalar type `T`.
"""
function initialize_buffers(::MonteCarlo, ::Type{T}, input::Engine.LineParameters, plan,
        buffers) where {T}
    dimensions = (size(observe(input, Engine.Z))..., plan.trials)
    Rs = Array{T}(undef, dimensions)
    Ls = similar(Rs)
    Cs = similar(Rs)
    Gs = similar(Rs)
    return merge(buffers, (; R = Rs, L = Ls, C = Cs, G = Gs))
end

function _record_sample!(
        buffers::NamedTuple{(:R, :L, :C, :G)},
        value::Engine.LineParameters,
        trial::Int,
        expected_frequencies
)
    impedance = observe(value, Engine.Z)
    size(impedance) == size(buffers.R)[1:3] || throw(DimensionMismatch(
        "Monte Carlo realizations produced incompatible impedance dimensions",
    ))
    observe(value, Engine.frequencies) == expected_frequencies || throw(DimensionMismatch(
        "Monte Carlo realizations produced incompatible frequency axes",
    ))
    for index in CartesianIndices(size(impedance))
        i, j, k = index.I
        buffers.R[i, j, k, trial] = observe(value, R, i, j, k)
        buffers.L[i, j, k, trial] = observe(value, L, i, j, k)
        buffers.G[i, j, k, trial] = observe(value, Engine.G, i, j, k)
        buffers.C[i, j, k, trial] = observe(value, C, i, j, k)
    end
    return buffers
end

function _map_samples(function_value, sample_values::Array{<:Real, 4})
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

function _aggregate(
        sample_values::NamedTuple{(:R, :L, :C, :G)},
        first_result::Engine.LineParameters,
        formulation::MonteCarlo
)
    summary_values = Tuple(
        _map_samples(SampleSummary, values)
    for values in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
    )
    summaries = NamedTuple{(:R, :L, :C, :G)}(summary_values)
    retained_details=_sample_shunt_model(first_result)
    source_details=details(first_result).data
    haskey(source_details,:comparison_unsupported) && (retained_details=ComputationDetails(
        merge(retained_details.data,(comparison_unsupported=source_details.comparison_unsupported,))))
    representation = materialize(first_result, summaries; details=retained_details)
    hist = formulation.options.data.return_histograms ?
           NamedTuple{(:R, :L, :C, :G)}(Tuple(
        _map_samples(
            values -> HistogramDensity(values; bins = formulation.options.data.bins),
            samples
        )
    for samples in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
    )) : nothing
    retained = formulation.options.data.return_samples ? sample_values : nothing
    return (; representation, statistics = summaries, samples = retained, histograms = hist)
end

_uncertain_arguments(point) = _uncertain_arguments(NamedTuple[], point, ())
# The uncertain arguments below `value` at `path`, appended to `records`.
function _uncertain_arguments(records, value::UncertainValue, path)
    return push!(records,(argument_path=path,nominal=value.nominal,standard_deviation=value.sigma))
end
function _uncertain_arguments(records, value::Gridpoint, path)
    for (index,argument) in enumerate(value.args)
        _uncertain_arguments(records,argument,(path...,index))
    end
    return records
end
function _uncertain_arguments(records, value::Union{Tuple,AbstractArray}, path)
    for (index,argument) in enumerate(value)
        _uncertain_arguments(records,argument,(path...,index))
    end
    return records
end
_uncertain_arguments(records, value, path) = records

# A draw is infeasible when its build throws `OverlapError`, or when its build or
# computation throws `DomainError`. Monte Carlo draws the trial again.
on_error(::MonteCarlo, exception::DomainError) = exception
on_error(::MonteCarlo, exception::DataModel.OverlapError) = exception

# Draw the arguments, build the problem and compute it with the point method of the inner
# formulation's backend: a trial is not a run. A draw whose exception `on_error` returns
# for the study is infeasible, and the next draw replaces it as if it was never drawn.
# `max_failures` infeasible draws in a row throw a `DomainError` with the uncertain
# arguments of the point: its distribution has no feasible draw. Any other exception
# propagates.
function _attempt(rng::Random.AbstractRNG, sampler)
    (; point, formulation) = sampler
    discarded = 0
    while true
        arguments = realize_arguments(rng, point, formulation.options.data.distribution)
        try
            return _attempt(realize(point, arguments), sampler)
        catch exception
            on_error(formulation, exception)
            discarded += 1
            discarded < formulation.options.data.max_failures || throw(DomainError(
                _uncertain_arguments(point),
                "Monte Carlo discarded $discarded infeasible draws in a row, the " *
                "max_failures limit; the last: " * sprint(showerror, exception)))
        end
    end
end

# The work of one trial. The type of a drawn problem is known only from its build, which is
# the user's: this method compiles the trial's computation from the problem it receives.
function _attempt(drawn::AbstractProblemDefinition, sampler)
    (; formulation, options) = sampler
    return compute(only(resolve(drawn, (formulation.inner,); options)), formulation.inner;
        options)
end

# A point of a Monte Carlo study builds its nominal problem from the means of its uncertain
# arguments before the first draw, and the constructor of the problem validates it.
materialize(::MonteCarlo, source) =
    realize(Random.Xoshiro(0), source, (_rng, mean, _sigma) -> mean)

# The run of a Monte Carlo study keeps the space of its problem: each point computes its
# trials from the uncertain point at its index.
traverse(space::ScalarOrGrid{AbstractProblemDefinition}, ::MonteCarlo, run::NamedTuple) =
    traverse(merge(run, (; space)))

# The input of a point of a Monte Carlo study is its nominal problem alone. The backend
# resolves each drawn problem of the point's trials.
resolve(::MonteCarlo, problem, formulations; options) =
    map(Returns(Resolved(problem, (;))), formulations)

# The value of a point of a Monte Carlo study comes from the trials of its uncertain point,
# the point of the problem's space at its index. The seed of the point combines the root
# seed of the study with that index. It returns the `MonteCarloResult` of that point, whose root
# seed is the seed of the point. `input.problem` is the nominal problem of the point. Each
# trial is a point of the inner formulation's backend, computed with the options of the
# problem that its owner normalized. A trial that ends in a later five-second interval of
# the point than the one in which it began logs the accepted count, the elapsed time and an
# estimate of the remaining time in the group `:progress`.
function traverse(study::MonteCarlo, run::NamedTuple, gridpoint::NamedTuple, input,
        ::AbstractFormulation)
    point = points(run.space, gridpoint.problem_index)
    seed = study.options.data.seed ⊻
           (UInt64(gridpoint.problem_index - 1) * 0x9e3779b97f4a7c15)
    rng = Random.Xoshiro(seed)
    progress = verbosity(run.record, :progress) > 0
    started = progress ? time_ns() : UInt64(0)
    retain_details = run.record.data.retain_details
    sampler = (; point, formulation = study, options = run.execution, retain_details)
    details_owner = typeof(study.inner)

    # Pilot: the first accepted realization.
    first_result = _attempt(rng, sampler)

    # The storage of every trial, allocated once from the pilot result.
    record = retain_details ? computation_details(details_owner, first_result) : nothing
    ntrials = something(
        study.options.data.trials,
        _dkw_trials(
            _observable_count(first_result),
            study.options.data.confidence,
            study.options.data.cdf_tol
        )
    )
    sample_values = initialize_buffers(study, Float64, first_result,
        (; trials = ntrials), (;))
    sample_axis = _sample_axis(first_result)
    retained = retain_details ? Vector{typeof(record)}(undef, ntrials) : nothing

    _trials!(rng, first_result, sample_values, retained, record, sampler, details_owner,
        sample_axis, ntrials, progress, started)
    aggregate = _aggregate(sample_values, first_result, study)
    # The trials computed with `study.inner`. The record of the first one holds that
    # formulation and the choices its point made, such as the equivalent earth of each earth
    # formula. The inputs are those of the nominal problem.
    completion = details(first_result).data[(:formulation, :choices)]
    physical_inputs = merge(Engine.completed_inputs(input.problem),
        (interpretation = :nominal_declaration, uncertain_arguments = _uncertain_arguments(point)))
    representation = Engine.retain_gridpoint(aggregate.representation, gridpoint;
        fields = merge(completion, (inputs = physical_inputs,
            uncertainty = (estimator = :empirical, representation = :marginal_mean_std),
            uncertainty_descriptions = _uncertainty_descriptions(MonteCarlo))))
    data = study.options.data
    return MonteCarloResult(data[(:confidence, :cdf_tol, :distribution)], [representation],
        fieldtype(ParametricResult, :failures)(), [aggregate.statistics],
        data.return_samples ? [aggregate.samples] : nothing,
        data.return_histograms ? [aggregate.histograms] : nothing, seed, [seed], [ntrials],
        retain_details ? ComputationDetails((; trials = [retained])) : ComputationDetails())
end

# The retained details of a point of a Monte Carlo study are those of its trials.
computation_details(::Type{<:AbstractFormulation}, point::MonteCarloResult) = details(point)

# Record the pilot result in `buffers`, then draw and record further trials until the study
# accepts `ntrials` of them. A trial that ends in a later five-second interval of the point
# than the one in which it began logs the accepted count, the elapsed time and an estimate
# of the remaining time: the remaining trials times the mean time of an accepted trial.
function _trials!(rng, first_result::R, buffers, retained, record, sampler, details_owner,
        sample_axis, ntrials::Int, progress::Bool, started::UInt64) where {R}
    accepted = 0
    value = first_result
    began = started
    while true
        accepted += 1
        _record_sample!(buffers, value, accepted, sample_axis)
        retained === nothing || (retained[accepted] = record)
        if progress
            now = time_ns()
            if fld(now - started, 5_000_000_000) > fld(began - started, 5_000_000_000)
                elapsed = (now - started) * 1e-9
                @info "Monte Carlo sampling progress" _group=:progress accepted trials=ntrials elapsed_seconds=elapsed eta_hours=(ntrials-accepted)*elapsed/accepted/3600
            end
            began = now
        end
        accepted < ntrials || break
        result = _attempt(rng, sampler)
        record = sampler.retain_details ? computation_details(details_owner, result) : nothing
        _sample_shunt_model(result) == _sample_shunt_model(first_result) || throw(ArgumentError(
            "Monte Carlo realizations changed shunt-model coverage; select a fixed geometry model for this study"))
        result isa R || throw(ArgumentError(
            "Monte Carlo realizations produced incompatible result types",
        ))
        value = result
    end
    return buffers
end

"""
    compute(problem::ParametricProblem, study::MonteCarlo; options)

Sample the uncertain arguments of each point of the problem's space and return the
[`MonteCarloResult`](@ref) of the computed points and the failed ones. The study runs
through [`traverse`](@ref LineCableModels.Commons.traverse): the problem's options are the
computation options of the inner formulation's backend for every trial, and `options` are
those of the run. It builds the nominal problem of every point before the first draw, then
computes the trials of each point from the point's seed: the root seed of the study, its
option `seed` or a fresh one that the computation fills in, combined with the position of
the point. The result joins the `MonteCarloResult` of each computed point in the order of
the space. A study with no computed point throws.
"""
function compute(problem::ParametricProblem, study::MonteCarlo;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt) === nothing &&
        throw(ArgumentError("MonteCarlo requires the Measurements extension to construct its result; " *
            "load it with `using Measurements` before compute"))
    # A study without a seed computes from a fresh root seed, which the study then stores.
    data = study.options.data
    seeded = if data.seed === nothing
        filled = merge(data, (seed = rand(Random.RandomDevice(), UInt64),))
        MonteCarlo(study.inner, FormulationOptions{typeof(filled)}(filled))
    else
        study
    end
    traversed = traverse(problem, seeded; options)
    # The result of each computed point, a `MonteCarloResult` of one point, joined in the
    # order of the space. The retained details of the run contain the details of the trials
    # of each point.
    computed = traversed.values
    retained = traversed.details.data
    return MonteCarloResult(
        data[(:confidence, :cdf_tol, :distribution)],
        [only(point.values) for point in computed],
        traversed.failures,
        [only(point.stats) for point in computed],
        data.return_samples ? [only(point.sample_values) for point in computed] : nothing,
        data.return_histograms ? [only(point.histogram_values) for point in computed] : nothing,
        seeded.options.data.seed,
        [only(point.point_seeds) for point in computed],
        [only(point.trial_counts) for point in computed],
        isempty(retained) ? ComputationDetails() :
        ComputationDetails((; trials = [only(point.data.trials) for point in retained.points]))
    )
end
