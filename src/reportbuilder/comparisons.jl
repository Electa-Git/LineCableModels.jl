"""
$(TYPEDSIGNATURES)

Select a completed scalar or formulation-space operand for comparison. Native
results supply their output terminal identities. An external result without that
information must be supplied as `(result=..., metadata=(port_order=..., ...))`.
"""
function select(definition::BenchmarkTableDefinition, result::AbstractCoreResult)
    coordinates=get(details(result), :coordinates, nothing)
    coordinates === nothing && throw(ArgumentError(
        "comparison needs output terminal identities; supply explicit operand metadata.port_order"))
    metadata=(port_order=coordinates,formulation=get(details(result),:formulations,(;)),axes=nothing)
    validate(definition,result,metadata)
    return (;result,metadata)
end

function select(definition::BenchmarkTableDefinition, result::ParametricResult)
    isempty(result.axes) && throw(ArgumentError("comparison requires problem/formulation axes"))
    first_operand=select(definition, first(result))
    for value in result
        select(definition, value).metadata.port_order == first_operand.metadata.port_order ||
            throw(ArgumentError("result-space output terminal identities differ"))
    end
    return (result, metadata=(port_order=first_operand.metadata.port_order,
        formulation=result.axes.formulations, axes=result.axes))
end

function select(definition::BenchmarkTableDefinition, operand::NamedTuple{(:result, :metadata)})
    ports=operand.metadata.port_order
    !isempty(ports) && allunique(ports) || throw(ArgumentError("output terminal identities must be nonempty and unique"))
    validate(definition,operand.result,operand.metadata)
    return operand
end

"""Select one retained result point without interpreting its formulation."""
function select(result::AbstractCoreResult, index::Integer)
    index == 1 || throw(ArgumentError("a scalar operand has only point 1"))
    return result
end
function select(result::ParametricResult, index::Integer)
    return result[index]
end

"""
$(TYPEDSIGNATURES)

Construct a comparison publication from explicitly ordered reference/candidate
operands. Raw operands are compared once. Supplying retained `comparisons` selects
those numerical products without recalculating them. `context` identifies the
case, benchmark and collection; it does not select physical equations.
"""
function select(definition::BenchmarkTableDefinition, source::NamedTuple)
    all(key -> haskey(source, key), (:reference, :candidate)) ||
        throw(ArgumentError("benchmark source requires an explicit reference and candidate"))
    isempty(setdiff(keys(source), (:reference, :candidate, :context, :comparisons, :settings))) ||
        throw(ArgumentError("unknown benchmark source fields"))
    reference=select(definition, source.reference)
    candidate=select(definition, source.candidate)
    reference.metadata.port_order == candidate.metadata.port_order ||
        throw(ArgumentError("reference and candidate output terminal identities differ"))
    context=get(source, :context, (id=:comparison, case_id=:standalone, collection=:manual))
    settings=get(source, :settings, definition.settings)
    if haskey(source, :comparisons)
        comparisons=source.comparisons
        isempty(comparisons) && throw(ArgumentError("benchmark has no retained comparisons"))
        for row in comparisons
            left=select(reference.result, row.reference_index)
            right=select(candidate.result, row.candidate_index)
            basis(left) === basis(right) && domain(left) === domain(right) ||
                throw(ArgumentError("retained comparison basis or domain differs"))
            size(row.error.absolute) == size(row.error.relative) ==
                (length(reference.metadata.port_order), length(reference.metadata.port_order)) ||
                throw(DimensionMismatch("RMS matrices differ from output terminal identities"))
        end
        isempty(definition.requested) && return (;reference,candidate,context,settings,comparisons)
        for key in definition.requested
            key in (:quantities,:statistics,:bands,:normalizations) && continue
            isequal(getproperty(definition.settings,key),getproperty(settings,key)) ||
                throw(ArgumentError("requested $key differs from the retained calculation; request explicit reanalysis"))
        end
        selected_settings=merge(settings,NamedTuple{Tuple(key for key in definition.requested if
            key in (:quantities,:statistics,:bands,:normalizations))}(Tuple(getproperty(definition.settings,key)
            for key in definition.requested if key in (:quantities,:statistics,:bands,:normalizations))))
        pairs=unique((row.reference_index,row.candidate_index) for row in comparisons)
        for quantity in selected_settings.quantities,statistic in selected_settings.statistics,
                band in selected_settings.bands,normalization in selected_settings.normalizations, pair in pairs
            any(row -> (row.reference_index,row.candidate_index)==pair && row.quantity==quantity &&
                row.statistic==statistic && row.error.details.band==band &&
                row.error.details.normalization==normalization,comparisons) ||
                throw(ArgumentError("requested comparison was not retained; request explicit reanalysis"))
        end
        comparisons=filter(row -> row.quantity in selected_settings.quantities &&
            row.statistic in selected_settings.statistics && row.error.details.band in selected_settings.bands &&
            row.error.details.normalization in selected_settings.normalizations,comparisons)
        return (; reference, candidate, context, settings=selected_settings, comparisons)
    end
    comparisons=select(definition,reference.result,candidate.result)
    return (; reference, candidate, context, settings, comparisons)
end

"""Calculate requested per-term comparisons through the result owner's compare methods."""
function select(definition::BenchmarkTableDefinition,
        reference::Union{AbstractCoreResult,ParametricResult},
        candidate::Union{AbstractCoreResult,ParametricResult})
    settings=definition.settings
    settings.statistics == (:value,) || throw(ArgumentError("line parameters require value comparisons"))
    if !(reference isa ParametricResult) && settings.pairing !== nothing
        count=candidate isa ParametricResult ? length(candidate) : 1
        length(settings.pairing)==count && all(pair -> first(pair)==1,settings.pairing) ||
            throw(ArgumentError("a scalar reference requires reference index 1 for every candidate"))
    end
    comparisons=NamedTuple[]
    for band in settings.bands, name in settings.quantities, normalization in settings.normalizations
        quantity=getproperty(Engine, name)
        controls=(; band, normalization, atol=settings.atol,
            fundamental=settings.fundamental, harmonics=settings.harmonics, unsupported=settings.unsupported)
        result=if reference isa ParametricResult
            candidate isa ParametricResult || throw(ArgumentError(
                "a reference space requires a candidate space and explicit pairing"))
            Engine.compare(reference, candidate, quantity; pairing=settings.pairing, controls...)
        else
            Engine.compare(reference, candidate, quantity; controls...)
        end
        errors=result isa ParametricResult ? result.values : (result,)
        for (index, error) in enumerate(errors)
            reference_index=settings.pairing === nothing ? 1 :
                first(only(filter(pair -> last(pair) == index, settings.pairing)))
            push!(comparisons, (quantity=name, statistic=:value, reference_index,
                candidate_index=index, error))
        end
    end
    return comparisons
end

"""
$(TYPEDSIGNATURES)

Expose comparison records, every matrix term and a labelled summary of per-term
maxima. Absolute RMS uses native units; relative RMS is displayed as percent.
Unavailable relative values retain their reasons and measured absolute values.
"""
function tabulate(definition::BenchmarkTableDefinition, source,
        published::NamedTuple{(:reference, :candidate, :context, :settings, :comparisons)})
    calculations=DataFrame([merge((;role), operand.metadata) for (role,operand) in
        pairs((reference=published.reference, candidate=published.candidate))])
    comparisons=NamedTuple[]
    terms=NamedTuple[]
    summaries=NamedTuple[]
    formulations=NamedTuple[]
    selections=map((published.reference,published.candidate)) do operand
        operand.result isa ParametricResult ?
            [value isa NamedTuple ? value : NamedTuple(value) for value in operand.result.axes.formulations] :
            [operand.metadata.formulation]
    end
    for (role, records) in ((:reference,selections[1]),(:candidate,selections[2]))
        labels=description(records;prefix=role === :reference ? "Reference F" : "F")
        for (index, record) in enumerate(records)
            push!(formulations, (role, formulation_index=index, label=labels[index], record))
        end
    end
    for (analysis,row) in enumerate(published.comparisons)
        error=row.error
        detail=error.details
        ports=published.candidate.metadata.port_order
        left=select(published.reference.result,row.reference_index)
        right=select(published.candidate.result,row.candidate_index)
        candidate_axes=published.candidate.result isa ParametricResult ? published.candidate.result.axes : nothing
        nproblems=candidate_axes === nothing ? 1 : length(candidate_axes.problems)
        problem_index=mod1(row.candidate_index,nproblems)
        formulation_index=cld(row.candidate_index,nproblems)
        unit=Units.native_unit(getproperty(Engine,row.quantity),basis(right))
        absolute=detach(error.absolute,1,definition.clip)
        relative=detach(error.relative,100,definition.clip)
        statuses=get(detail,:status,fill(:not_recorded,size(absolute)))
        local_reasons=get(detail,:normalization_reason,fill(nothing,size(absolute)))
        reasons=map(reason -> reason === nothing ? get(detail,:reason,nothing) : reason,local_reasons)
        identity=(analysis, snapshot=get(row,:snapshot,"live"), benchmark=string(published.context.id), case_id=string(published.context.case_id),
            collection=string(published.context.collection), quantity=row.quantity, statistic=row.statistic,
            reference_point=row.reference_index, candidate_point=row.candidate_index,
            problem_index, formulation_index, band=detail.band, normalization=detail.normalization,
            absolute_unit=Units.label(unit),samples=detail.sample_count,
            requested_bounds_Hz=get(detail,:requested_bounds,missing), actual_bounds_Hz=detail.actual_bounds)
        push!(comparisons,merge(identity,(absolute_rms=absolute,relative_rms_percent=relative,
            status=copy(statuses),reason=reasons,port_order=copy(ports),
            sample_indices=copy(get(detail,:indices,Int[])),tolerance=get(detail,:atol,missing))))
        for i in eachindex(ports),j in eachindex(ports)
            push!(terms,merge(identity,(row=i,column=j,response=ports[i],excitation=ports[j],
                absolute_rms=absolute[i,j],relative_rms_percent=relative[i,j],
                status=statuses[i,j],reason=reasons[i,j])))
        end
        ai=findall(!ismissing,absolute)
        ri=findall(!ismissing,relative)
        amax=isempty(ai) ? nothing : ai[argmax(absolute[ai])]
        rmax=isempty(ri) ? nothing : ri[argmax(relative[ri])]
        push!(summaries,merge(identity,(
            maximum_absolute_rms=amax === nothing ? missing : absolute[amax],
            absolute_term=amax === nothing ? missing : (ports[amax[1]],ports[amax[2]]),
            maximum_relative_rms_percent=rmax === nothing ? missing : relative[rmax],
            relative_term=rmax === nothing ? missing : (ports[rmax[1]],ports[rmax[2]]),
            unavailable=count(ismissing,relative), term_count=length(relative),
            reasons=unique(filter(!isnothing,vec(reasons))))))
    end
    return (calculations, formulations=DataFrame(formulations),comparisons=DataFrame(comparisons),
        terms=DataFrame(terms),maxima=DataFrame(summaries),summary=tabulate(definition,source,summaries))
end

function illustrate(definition::BenchmarkTableDefinition, source,
        published::NamedTuple{(:reference, :candidate, :context, :settings, :comparisons)}, table)
    illustration=definition.illustration
    (illustration === nothing || illustration === false) && return nothing
    illustration === true && return PlotBuilder.plot(published; definition.plot_options...)
    return illustration(published; definition.plot_options...)
end

"""
$(TYPEDSIGNATURES)

Format retained per-term maxima as one compact row per case, problem, formulation,
band and normalization. Quantity columns retain unavailable counts and identify
the largest term; no averaging across terms, cases or formulations is performed.
"""
function tabulate(::BenchmarkTableDefinition,source,maxima::AbstractVector{<:NamedTuple})
    rows=Dict{Tuple,Dict{Symbol,Any}}()
    order=Tuple[]
    for metric in maxima
        key=(get(metric,:snapshot,"legacy"),metric.case_id,metric.benchmark,metric.collection,metric.problem_index,
            metric.formulation_index,metric.reference_point,metric.band,metric.normalization,metric.statistic)
        if !haskey(rows,key)
            push!(order,key)
            rows[key]=Dict(:snapshot=>get(metric,:snapshot,"legacy"),:case_id=>metric.case_id,:benchmark=>metric.benchmark,:collection=>metric.collection,
                :problem_index=>metric.problem_index,:formulation_index=>metric.formulation_index,
                :reference_point=>metric.reference_point,:band=>metric.band,:normalization=>metric.normalization,
                :statistic=>metric.statistic,:samples=>metric.samples,
                :requested_bounds_Hz=>metric.requested_bounds_Hz,:actual_bounds_Hz=>metric.actual_bounds_Hz)
        end
        text=if !ismissing(metric.maximum_relative_rms_percent)
            string(round(metric.maximum_relative_rms_percent;sigdigits=5),"% · ",metric.relative_term)
        elseif !ismissing(metric.maximum_absolute_rms)
            string("relative unavailable; absolute ",round(metric.maximum_absolute_rms;sigdigits=5),
                " ",metric.absolute_unit," · ",metric.absolute_term)
        else
            "unavailable (" * join(string.(metric.reasons),"; ") * ")"
        end
        metric.unavailable > 0 && (text *= "; $(metric.unavailable)/$(metric.term_count) unavailable")
        rows[key][metric.quantity]=text
    end
    frame=DataFrame()
    for key in order
        push!(frame,rows[key];cols=:union)
    end
    return frame
end

"""Describe a report frequency band in ordinary scientific terms."""
function description(::BenchmarkTableDefinition,band)
    names=(all="Entire range",dc="Near DC",harmonic="Harmonic range",narrow="Narrowband",wide="Wideband")
    return band isa Symbol ? get(names,band,string(band)) : string(band," Hz")
end

"""Check the explicit output coordinates before a report accesses matrix entries."""
function validate(::BenchmarkTableDefinition,result::AbstractCoreResult,metadata::NamedTuple)
    ports=metadata.port_order
    size(observe(result,Z))[1:2] == (length(ports),length(ports)) ||
        throw(DimensionMismatch("output terminal identities do not match the matrix dimensions"))
    get(details(result),:coordinates,ports) == ports || throw(ArgumentError("explicit output terminal identities differ from the result"))
    return nothing
end
function validate(definition::BenchmarkTableDefinition,result::ParametricResult,metadata::NamedTuple)
    isempty(result.axes) && throw(ArgumentError("result-space axes are required"))
    foreach(value -> validate(definition,value,metadata),result)
    return nothing
end
