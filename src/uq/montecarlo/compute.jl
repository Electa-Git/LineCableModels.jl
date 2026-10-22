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

_uncertain_arguments(point::Gridpoint) = _uncertain_arguments(NamedTuple[], point, ())
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

"""
    compute(resolved::Resolved, formulation::MonteCarlo; options)

Compute the trials of one point of a Monte Carlo study and return the `MonteCarloResult` of
that point, whose root seed is the seed of the point. `resolved.problem` is the nominal
problem of the point, built from the means of its uncertain arguments before the first
draw, and `resolved.input` holds the uncertain point, its seed and its gridpoint identity.
Each trial is a point of the inner formulation's backend, computed with the `options` that
its owner normalized. A trial that ends in a later five-second interval of the point than
the one in which it began logs the accepted count, the elapsed time and an estimate of the
remaining time in the group `:progress`.
"""
function compute(resolved::Resolved, formulation::MonteCarlo; options::ComputationOptions)
    (; point, seed, gridpoint) = resolved.input
    rng = Random.Xoshiro(seed)
    progress = verbosity(options, :progress) > 0
    started = progress ? time_ns() : UInt64(0)
    sampler = (; point, formulation, options)
    details_owner = typeof(formulation.inner)

    # Pilot: the first accepted realization.
    first_result = _attempt(rng, sampler)

    # The storage of every trial, allocated once from the pilot result.
    record = formulation.options.data.retain_details ?
             computation_details(details_owner, first_result) : nothing
    ntrials = something(
        formulation.options.data.trials,
        _dkw_trials(
            _observable_count(first_result),
            formulation.options.data.confidence,
            formulation.options.data.cdf_tol
        )
    )
    sample_values = initialize_buffers(formulation, Float64, first_result,
        (; trials = ntrials), (;))
    sample_axis = _sample_axis(first_result)
    retained = formulation.options.data.retain_details ?
               Vector{typeof(record)}(undef, ntrials) : nothing

    _trials!(rng, first_result, sample_values, retained, record, sampler, details_owner,
        sample_axis, ntrials, progress, started)
    aggregate = _aggregate(sample_values, first_result, formulation)
    # The trials computed with `formulation.inner`. The record of the first one holds that
    # formulation and the choices its point made, such as the equivalent earth of each earth
    # formula. The inputs are those of the nominal problem.
    completion = details(first_result).data[(:formulation, :choices)]
    physical_inputs = merge(Engine.completed_inputs(resolved.problem),
        (interpretation = :nominal_declaration, uncertain_arguments = _uncertain_arguments(point)))
    representation = Engine.retain_gridpoint(aggregate.representation, gridpoint;
        fields = merge(completion, (inputs = physical_inputs,
            uncertainty = (estimator = :empirical, representation = :marginal_mean_std),
            uncertainty_descriptions = _uncertainty_descriptions(MonteCarlo))))
    data = formulation.options.data
    return MonteCarloResult(data[(:confidence, :cdf_tol, :distribution)], [representation],
        fieldtype(ParametricResult, :failures)(), [aggregate.statistics],
        data.return_samples ? [aggregate.samples] : nothing,
        data.return_histograms ? [aggregate.histograms] : nothing, seed, [seed], [ntrials],
        data.retain_details ? ComputationDetails((; trials = [retained])) : ComputationDetails())
end

# Record the pilot result in `buffers`, then draw and record further trials until the study
# accepts `ntrials` of them. A trial that ends in a later five-second interval of the point
# than the one in which it began logs the accepted count, the elapsed time and an estimate
# of the remaining time: the remaining trials times the mean time of an accepted trial.
function _trials!(rng, first_result::R, buffers, retained, record, sampler, details_owner,
        sample_axis, ntrials::Int, progress::Bool, started::UInt64) where {R}
    formulation = sampler.formulation
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
        record = formulation.options.data.retain_details ?
                 computation_details(details_owner, result) : nothing
        _sample_shunt_model(result) == _sample_shunt_model(first_result) || throw(ArgumentError(
            "Monte Carlo realizations changed shunt-model coverage; select a fixed geometry model for this study"))
        result isa R || throw(ArgumentError(
            "Monte Carlo realizations produced incompatible result types",
        ))
        value = result
    end
    return buffers
end

# A Monte Carlo study builds the nominal problem of each point of its space from the means
# of the uncertain arguments, before the first draw: the constructor of the problem
# validates it. When `on_error` returns the exception of a nominal build for the
# `Combinatorial` study of the points, the point fails and the run continues. Each built
# point is resolved here, with the uncertain point, its seed and its gridpoint identity,
# and the runner computes its trials as one point.
function traverse(problem::ScalarOrGrid{AbstractProblemDefinition}, formulation::MonteCarlo,
        root_seed::UInt64; options::ComputationOptions)
    sources = points(problem)
    formulations = [formulation]
    # The owner of the inner formulation normalizes the options of the trials once.
    execution = computation_options(typeof(formulation.inner), options)
    levels = verbosity(get(options.data, :verbosity, (default = 0,)))
    on_result = get(options.data, :on_result, nothing)
    gridpoints = CartesianIndices((length(sources), 1))
    source_id = Commons.gridpoint_id().source_id
    progress = get(levels, :progress, levels.default) > 0
    timing = get(levels, :timing, levels.default) > 0
    started = time_ns()
    scoped = haskey(options.data, :verbosity) || on_result !== nothing
    study = Combinatorial(formulation)
    failures = fieldtype(ParametricResult, :failures)()
    attempts = map(enumerate(sources)) do (problem_index, source)
        gridpoint = Commons.gridpoint_id(; source_id, problem_index)
        try
            nominal = realize(Random.Xoshiro(0), source, (_rng, mean, _sigma) -> mean)
            Resolved(nominal, (; point = source,
                seed = root_seed ⊻ (UInt64(problem_index - 1) * 0x9e3779b97f4a7c15),
                gridpoint))
        catch exception
            on_error(study, exception)
            push!(failures, gridpoint => exception)
            exception
        end
    end
    problem_indices = isempty(failures) ? eachindex(attempts) :
                      setdiff!(collect(eachindex(attempts)),
        (gridpoint.problem_index for (gridpoint, _) in failures))
    problems = isempty(failures) ? attempts :
               map(Base.Fix1(getindex, attempts), problem_indices)
    return Base.invokelatest(traverse, problem, formulation,
        (; problems, problem_indices, study, failures, formulations, execution, on_result,
            gridpoints, spaced = true, source_id, progress, timing, started), scoped, levels)
end

# Monte Carlo samples the arguments of the points of a space. A scalar problem is the space
# of one point, whose argument is the problem.
function compute(problem::ParametricProblem{<:AbstractProblemDefinition}, formulation::MonteCarlo)
    space = Gridspace{typeof(problem.space)}(identity, (Grid((problem.space,)),))
    return compute(ParametricProblem(space, problem.options), formulation)
end

"""
    compute(problem::ParametricProblem, formulation::MonteCarlo)

Sample the uncertain arguments of each point of the problem's space and return the
[`MonteCarloResult`](@ref) of the computed points and the failed ones. The study runs
through [`traverse`](@ref LineCableModels.Commons.traverse) with the problem's options: it
builds the nominal problem of every point before the first draw, then computes each point
with `compute(resolved, formulation; options)`, from the point's seed, the root seed
`formulation.options.data.seed` or a fresh one, combined with the position of the point.
The result joins the `MonteCarloResult` of each computed point in the order of the space.
A study with no computed point throws.
"""
function compute(problem::ParametricProblem, formulation::MonteCarlo)
    Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt) === nothing &&
        throw(ArgumentError("MonteCarlo requires the Measurements extension to construct its result; " *
            "load it with `using Measurements` before compute"))
    data = formulation.options.data
    root_seed = data.seed === nothing ? rand(Random.RandomDevice(), UInt64) : data.seed
    traversed = traverse(problem.space, formulation, root_seed; options = problem.options)
    # The result of each computed point, a `MonteCarloResult` of one point, joined in the
    # order of the space.
    computed = traversed.values
    return MonteCarloResult(
        data[(:confidence, :cdf_tol, :distribution)],
        [only(point.values) for point in computed],
        traversed.failures,
        [only(point.stats) for point in computed],
        data.return_samples ? [only(point.sample_values) for point in computed] : nothing,
        data.return_histograms ? [only(point.histogram_values) for point in computed] : nothing,
        root_seed,
        [only(point.point_seeds) for point in computed],
        [only(point.trial_counts) for point in computed],
        data.retain_details ? ComputationDetails((;
            trials = [only(details(point).data.trials) for point in computed])) :
        ComputationDetails()
    )
end
