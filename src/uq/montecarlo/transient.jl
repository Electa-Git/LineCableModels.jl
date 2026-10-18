function _sample_shunt_model(value::Engine.TransientAnalysis.TransientResponse)
    line = get(details(value).data,:line,ComputationDetails())
    return _sample_shunt_model(line)
end

function _observable_count(value::Engine.TransientAnalysis.TransientResponse)
    return length(value.voltage) + length(value.current) + 8prod(size(value.voltage)[1:2])
end

function _sample_axis(value::Engine.TransientAnalysis.TransientResponse, point, formulation, child_options)
    return (; time=value.t, coordinates=value.coordinates)
end

function initialize_buffers(::Union{MonteCarlo,LinearError}, ::Type, input::Engine.TransientAnalysis.TransientResponse, plan, buffers)
    T = eltype(input.voltage)
    samples = (voltage=Array{T}(undef,size(input.voltage)...,plan.trials),
        current=Array{T}(undef,size(input.current)...,plan.trials))
    kpis = map((:voltage,:current)) do _
        (peak=Array{T}(undef,size(input.voltage)[1:2]...,plan.trials),
            peak_time=Array{T}(undef,size(input.voltage)[1:2]...,plan.trials),
            arrival_time=Array{T}(undef,size(input.voltage)[1:2]...,plan.trials),
            final_value=Array{T}(undef,size(input.voltage)[1:2]...,plan.trials))
    end
    return merge(buffers,samples,(kpis=NamedTuple{(:voltage,:current)}(kpis),))
end

function _record_sample!(buffers::NamedTuple, value::Engine.TransientAnalysis.TransientResponse, trial::Int, axis)
    value.t == axis.time && value.coordinates == axis.coordinates ||
        throw(DimensionMismatch("Monte Carlo transient realizations must share time and terminal axes"))
    size(value.voltage) == size(buffers.voltage)[1:3] ||
        throw(DimensionMismatch("Monte Carlo transient realizations produced incompatible response dimensions"))
    buffers.voltage[:,:,:,trial] .= value.voltage
    buffers.current[:,:,:,trial] .= value.current
    for (name,selector) in ((:voltage,Engine.TransientAnalysis.voltage),(:current,Engine.TransientAnalysis.current))
        peak = Engine.TransientAnalysis.peak(value,selector)
        destination = getproperty(buffers.kpis,name)
        destination.peak[:,:,trial] .= peak.value
        destination.peak_time[:,:,trial] .= peak.time
        destination.arrival_time[:,:,trial] .= Engine.TransientAnalysis.arrival_time(value,selector)
        destination.final_value[:,:,trial] .= Engine.TransientAnalysis.dc_offset(value,selector)
    end
    return buffers
end

function _aggregate(buffers::NamedTuple, source::Engine.TransientAnalysis.TransientResponse,
        formulation::MonteCarlo)
    # An identically zero trace has no arrival time. Preserve that missing KPI
    # rather than inventing an arrival or dropping nonzero terminal samples.
    summarize = values -> all(isfinite,values) ? SampleSummary(values) : missing
    summaries = (voltage=_map_samples(SampleSummary,buffers.voltage),
        current=_map_samples(SampleSummary,buffers.current),
        kpis=map(quantities -> map(values -> _map_samples(summarize,values),quantities),buffers.kpis))
    representation = materialize(source, summaries, formulation)
    histograms = if formulation.options.data.return_histograms
        density = values -> all(isfinite,values) ?
            HistogramDensity(values; bins=formulation.options.data.bins) : missing
        (voltage=_map_samples(density,buffers.voltage),current=_map_samples(density,buffers.current),
            kpis=map(quantities -> map(values -> _map_samples(density,values),quantities),buffers.kpis))
    else
        nothing
    end
    return (; representation, statistics=summaries,
        samples=formulation.options.data.return_samples ? buffers : nothing, histograms)
end

function materialize(source::Engine.TransientAnalysis.TransientResponse,
        quantities::NamedTuple{(:voltage,:current,:kpis)}, ::Union{LinearError,MonteCarlo})
    voltages = materialize.(quantities.voltage)
    currents = materialize.(quantities.current)
    T = eltype(voltages)
    return Engine.TransientAnalysis.TransientResponse(convert(Vector{T},source.t),voltages,currents,
        source.coordinates,details(source))
end

function validate(value::Engine.TransientAnalysis.TransientResponse, ::Type{MonteCarloResult},
        statistics_product, sample_product, histogram_product, trials)
    expected = (:voltage,:current,:kpis)
    validate(statistics_product,MonteCarloResult,expected,"transient statistics products")
    sample_product === nothing || validate(sample_product,MonteCarloResult,expected,"transient sample products")
    histogram_product === nothing || validate(histogram_product,MonteCarloResult,expected,"transient histogram products")
    for name in (:voltage,:current)
        quantity = getproperty(value,name)
        summaries = getproperty(statistics_product,name)
        Engine.has_uncertainty_type(eltype(quantity)) || throw(ArgumentError(
            "MonteCarloResult requires uncertainty-bearing transient samples"))
        size(summaries) == size(quantity) && all(summary -> summary.n == trials,summaries) ||
            throw(DimensionMismatch("transient summaries must match end, terminal, time and trial axes"))
        sample_product === nothing || size(getproperty(sample_product,name)) == (size(quantity)...,trials) ||
            throw(DimensionMismatch("transient samples must match end, terminal, time and trial axes"))
    end
    return value
end
