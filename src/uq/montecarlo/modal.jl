# A sampled line and its modal decomposition form one trial, before aggregation.
struct ModalComposition{L,M,O} <: AbstractFormulation
    line::L
    modal::M
    options::O
end

function compute(problem::Engine.LineParametersProblem, formulation::ModalComposition;
        options::Union{NamedTuple,ComputationOptions}=ComputationOptions())
    phase = compute(problem, formulation.line; options)
    return compute(LineCableModels.ModalAnalysisProblem(phase), formulation.modal;
        options=formulation.options)
end

function compute(problem::ParametricProblem, formulation::MonteCarlo,
        modal::LineCableModels.ModalAnalysisFormulation;
        modal_options::Union{NamedTuple,ComputationOptions}=ComputationOptions())
    composed = ModalComposition(formulation.inner, modal, ComputationOptions(modal_options))
    sampling = MonteCarlo(composed, ComputationOptions(merge(formulation.options.data,
        (retain_details=true,))))
    return compute(problem, sampling)
end

function Base.NamedTuple(formulation::ModalComposition)
    return (line=NamedTuple(formulation.line),modal=NamedTuple(formulation.modal),
        options=formulation.options.data)
end

function computation_details(::Type{<:ModalComposition}, result::Engine.LineParameters)
    return details(result)
end

function _observable_count(value::Engine.LineParameters{<:Complex,<:Real,<:Engine.ModalDomain})
    n, _, nf = size(value.Z)
    return (8n*n + 9n)*nf
end

const _MODAL_OVERLAP_THRESHOLD = 0.8
struct ModeAssignmentError <: Exception
    overlap::Matrix{Float64}
end
function Base.showerror(io::IO, error::ModeAssignmentError)
    print(io, "modal current-basis overlap ", error.overlap,
        " falls below the mode-assignment threshold ", _MODAL_OVERLAP_THRESHOLD,
        "; this trial was excluded without replacement")
end

function _sample_axis(value::Engine.LineParameters{<:Complex,<:Real,<:Engine.ModalDomain}, point,
        formulation::MonteCarlo{<:ModalComposition}, child_options)
    declaration = realize(Random.Xoshiro(0), point, (_rng, mean, _sigma) -> mean)
    reference = compute(declaration, formulation.inner; options=child_options)
    return (; frequencies=frequencies(reference), reference)
end

function initialize_buffers(::Union{MonteCarlo,LinearError}, ::Type{T},
        input::Engine.LineParameters{<:Complex,<:Real,<:Engine.ModalDomain}, plan, buffers) where {T}
    segment = LineCableModels.PropagationParameters(input)
    quantities = (R=R(input), L=L(input), C=C(input), G=Engine.G(input),
        gamma=LineCableModels.gamma(input), Zc=LineCableModels.Zc(input),
        Yc=LineCableModels.Yc(input), Ti=LineCableModels.Ti(input), Tv=LineCableModels.Tv(input),
        alpha=LineCableModels.alpha(segment),
        beta=LineCableModels.beta(segment),
        velocity=LineCableModels.velocity(segment))
    products = map(values(quantities)) do quantity
        Array{eltype(quantity)}(undef, size(quantity)..., plan.trials)
    end
    return merge(buffers, (products=NamedTuple{keys(quantities)}(products),
        alignment=initialize_buffers(Engine.ModalAnalysis.align_modes,T,input,(;),(;)),
        accepted=Ref(0),
        cluster_sizes=Matrix{Int}(undef,size(input.Z,1),size(input.Z,3)),
        overlaps=Array{T}(undef,size(input.Z,1),size(input.Z,3),plan.trials),
        permutations=Array{Int}(undef,size(input.Z,1),size(input.Z,3),plan.trials)))
end

function _record_sample!(buffers::NamedTuple, value::Engine.LineParameters{<:Complex,<:Real,<:Engine.ModalDomain},
        trial::Int, axis)
    assignment = Engine.ModalAnalysis.align_modes(value,axis.reference;
        buffers=buffers.alignment)
    buffers.overlaps[:,:,trial] .= assignment.overlap
    buffers.permutations[:,:,trial] .= assignment.permutation
    buffers.cluster_sizes .= assignment.cluster_sizes
    minimum(assignment.overlap) >= _MODAL_OVERLAP_THRESHOLD ||
        return false
    value = assignment.parameters
    segment = LineCableModels.PropagationParameters(value)
    quantities = (R=R(value), L=L(value), C=C(value), G=Engine.G(value),
        gamma=LineCableModels.gamma(value), Zc=LineCableModels.Zc(value),
        Yc=LineCableModels.Yc(value), Ti=LineCableModels.Ti(value), Tv=LineCableModels.Tv(value),
        alpha=LineCableModels.alpha(segment), beta=LineCableModels.beta(segment),
        velocity=LineCableModels.velocity(segment))
    for name in keys(quantities)
        quantity = getproperty(quantities,name)
        destination = selectdim(getproperty(buffers.products,name),ndims(quantity)+1,trial)
        destination .= quantity
    end
    buffers.accepted[] = trial
    return true
end

function _map_samples(function_value, sample_values::AbstractArray{<:Complex})
    return (real=_map_samples(function_value, real.(sample_values)),
        imag=_map_samples(function_value, imag.(sample_values)))
end

function _aggregate(buffers::NamedTuple, source::Engine.LineParameters{<:Complex,<:Real,<:Engine.ModalDomain},
        formulation::MonteCarlo)
    trials = buffers.accepted[]
    products = map(values(buffers.products)) do quantity
        copy(selectdim(quantity,ndims(quantity),1:trials))
    end
    samples = NamedTuple{keys(buffers.products)}(products)
    summaries = map(quantity -> _map_samples(SampleSummary,quantity),samples)
    histograms = formulation.options.data.return_histograms ? map(samples) do quantity
        _map_samples(values -> HistogramDensity(values; bins=formulation.options.data.bins),quantity)
    end : nothing
    quantities = merge(summaries, map((; gamma=summaries.gamma, Tv=summaries.Tv,
        Ti=summaries.Ti)) do quantity
        complex.(materialize.(quantity.real),materialize.(quantity.imag))
    end)
    representation = materialize(source, quantities, formulation)
    assignment = (threshold=_MODAL_OVERLAP_THRESHOLD, frequencies=copy(frequencies(source)),
        cluster_sizes=copy(buffers.cluster_sizes),
        overlap=copy(buffers.overlaps[:,:,1:trials]),
        permutation=copy(buffers.permutations[:,:,1:trials]))
    return (; representation, statistics=summaries,
        samples=formulation.options.data.return_samples ? samples : nothing,
        histograms, mode_assignment=assignment)
end

# Both propagation methods construct the modal core from these physical products.
function materialize(source::Engine.LineParameters{<:Complex,<:Real,<:Engine.ModalDomain},
        quantities::NamedTuple{(:R,:L,:C,:G,:gamma,:Zc,:Yc,:Ti,:Tv,:alpha,:beta,:velocity)},
        ::Union{LinearError,MonteCarlo})
    domain = Engine.ModalDomain(LineCableModels.ModalOperators(
        quantities.Tv, quantities.Ti), quantities.gamma)
    angular = reshape(2π .* frequencies(source),1,1,:)
    return Engine.LineParameters(domain,
        complex.(materialize.(quantities.R),materialize.(quantities.L).*angular),
        complex.(materialize.(quantities.G),materialize.(quantities.C).*angular),
        frequencies(source); basis=basis(source), details=details(source))
end

function validate(value::Engine.LineParameters{<:Complex,<:Real,<:Engine.ModalDomain}, ::Type{MonteCarloResult},
        statistics_product, sample_product, histogram_product, trials)
    expected = (:R,:L,:C,:G,:gamma,:Zc,:Yc,:Ti,:Tv,:alpha,:beta,:velocity)
    validate(statistics_product,MonteCarloResult,expected,"modal statistics products")
    sample_product === nothing || validate(sample_product,MonteCarloResult,expected,"modal sample products")
    histogram_product === nothing || validate(histogram_product,MonteCarloResult,expected,"modal histogram products")
    Engine.has_uncertainty_type(eltype(value.Z)) ||
        throw(ArgumentError("MonteCarloResult requires uncertainty-bearing modal coefficients"))
    n, _, nf = size(value.Z)
    for name in expected
        shape = name in (:R,:L,:C,:G,:Ti,:Tv) ? (n,n,nf) : (n,nf)
        product = getproperty(statistics_product,name)
        components = name in (:gamma,:Zc,:Yc,:Ti,:Tv) ? values(product) : (product,)
        for component in components
            size(component) == shape && all(summary -> summary.n == trials,component) ||
                throw(DimensionMismatch("modal summaries must match mode, frequency and trial axes"))
        end
        sample_product === nothing || size(getproperty(sample_product,name)) == (shape...,trials) ||
            throw(DimensionMismatch("modal samples must match mode, frequency and trial axes"))
    end
    return value
end

# A failed modal correspondence consumes its trial and is recorded without replacement.
function _trials!(rng, first_result::R, buffers, retained, record, sampler, details_owner,
        sample_axis, trial_timings, attempts::Int, ntrials::Int, progress::Bool,
        started::UInt64) where {R<:Engine.LineParameters{<:Complex,<:Real,<:Engine.ModalDomain}}
    formulation = sampler.formulation
    failures = sampler.failures
    previous = started
    last_log = started
    average_seconds = 0.0
    accepted = 0
    unmatched = 0
    value = first_result
    while true
        matched = _record_sample!(buffers, value, accepted + 1, sample_axis)
        if !matched
            exception = ModeAssignmentError(copy(@view buffers.overlaps[:,:,accepted+1]))
            push!(failures, _failure_record(attempts, accepted + 1, :compute,
                sampler.drawn[], exception, backtrace()))
            unmatched += 1
        end
        if matched
            accepted += 1
            retained === nothing || (retained[accepted] = record)
            trial_timings === nothing || push!(trial_timings, details(value).data.timing)
        end
        if progress
            now = time_ns()
            interval = (now - previous) * 1e-9
            average_seconds = accepted == 1 ? interval : 0.2 * interval + 0.8 * average_seconds
            previous = now
            if now - last_log >= 5_000_000_000
                @info "Monte Carlo sampling progress" _group=:progress accepted trials=ntrials attempts rejected=length(failures) elapsed_seconds=(now-started)*1e-9 eta_hours=(ntrials-accepted)*average_seconds/3600
                last_log = now
            end
        end
        accepted + unmatched < ntrials || break
        result = nothing
        while result === nothing
            attempts += 1
            result = _attempt(rng, sampler, attempts, accepted)
        end
        record = formulation.options.data.retain_details ?
                 computation_details(details_owner, result) : nothing
        _sample_shunt_model(result) == _sample_shunt_model(first_result) || throw(ArgumentError(
            "Monte Carlo realizations changed shunt-model coverage; select a fixed geometry model for this study"))
        typeof(result) === R || throw(ArgumentError(
            "Monte Carlo realizations produced incompatible result types",
        ))
        if retained !== nothing
            typeof(record) === eltype(retained) || throw(ArgumentError(
                "Monte Carlo trials produced incompatible details record types",
            ))
        end
        value = result
    end
    return (; attempts, accepted, last_log)
end
