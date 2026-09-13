"""
$(TYPEDSIGNATURES)

Select an owned result and its explicit output identities. External results
without terminal metadata must be supplied as `(result, metadata)`.
"""
function select(definition::BenchmarkTableDefinition, result::AbstractCoreResult)
    coordinates=get(details(result), :coordinates, nothing)
    coordinates === nothing && throw(ArgumentError(
        "comparison needs output terminal identities; supply explicit operand metadata.port_order"))
    metadata=(port_order=coordinates,formulation=get(details(result),:formulations,(;)),axes=nothing)
    validate(definition,result,metadata)
    return select(definition,(;result,metadata))
end

function select(definition::BenchmarkTableDefinition, result::ParametricResult)
    isempty(result.axes) && throw(ArgumentError("comparison requires problem/formulation axes"))
    first_operand=select(definition, first(result))
    metadata=(port_order=first_operand.metadata.port_order,
        formulation=result.axes.formulations, axes=result.axes)
    validate(definition,result,metadata)
    return select(definition,(;result,metadata))
end

function select(definition::BenchmarkTableDefinition, result::AbstractUncertaintyResult)
    isempty(result) && throw(ArgumentError("cannot compare an empty uncertainty result"))
    metadata=merge(select(definition, first(result)).metadata,
        (formulation=NamedTuple(result).formulation,))
    validate(definition,result,metadata)
    return select(definition,(;result,metadata))
end

function select(definition::BenchmarkTableDefinition, operand::NamedTuple{(:result, :metadata)})
    ports=operand.metadata.port_order
    !isempty(ports) && allunique(ports) || throw(ArgumentError("output terminal identities must be nonempty and unique"))
    validate(definition,operand.result,operand.metadata)
    raw=operand.result isa ParametricResult ? operand.result.axes.formulations :
        operand.result isa AbstractUncertaintyResult ? [NamedTuple(operand.result).formulation] : [operand.metadata.formulation]
    sources=map(raw) do value
        LineCableModels.ImportExport.deserialize_value(Val(:formulation),value)
    end
    return (result=operand.result,metadata=merge(operand.metadata,(formulation_sources=sources,)))
end

"""Select a core point without selecting a statistic or pooling a population."""
function select(result::AbstractCoreResult, index::Integer)
    index == 1 || throw(ArgumentError("a scalar operand has only point 1"))
    return result
end
select(result::Union{ParametricResult,AbstractUncertaintyResult}, index::Integer) = result[index]
function select(result::ObservationPublication,index::Integer)
    index==1 || throw(ArgumentError("a detached publication retains one selected point"))
    return result
end

"""
$(TYPEDSIGNATURES)

Select explicit operands and either compare them once or render retained
comparisons unchanged. Measurement evidence is supplied as completed records;
this action does not open evidence files or run calculations.
"""
function select(definition::BenchmarkTableDefinition, source::NamedTuple)
    all(key -> haskey(source,key),(:reference,:candidate)) ||
        throw(ArgumentError("benchmark source requires an explicit reference and candidate"))
    isempty(setdiff(keys(source),(:reference,:candidate,:context,:comparisons,:settings,:measurements))) ||
        throw(ArgumentError("unknown benchmark source fields"))
    reference=select(definition,source.reference)
    candidate=select(definition,source.candidate)
    reference.metadata.port_order == candidate.metadata.port_order ||
        throw(ArgumentError("reference and candidate output terminal identities differ"))
    context=get(source,:context,(id=:comparison,case_id=:standalone,collection=:manual))
    measurements=get(source,:measurements,nothing)
    settings=get(source,:settings,definition.settings)
    settings=BenchmarkTableDefinition(;settings...).settings
    if haskey(source,:comparisons)
        comparisons=source.comparisons
        isempty(comparisons) && throw(ArgumentError("benchmark has no retained comparisons"))
        for row in comparisons
            left=select(reference.result,row.reference_index)
            right=select(candidate.result,row.candidate_index)
            left_domain=left isa ObservationPublication ? reference.metadata.domain : nameof(domain(left))
            right_domain=right isa ObservationPublication ? candidate.metadata.domain : nameof(domain(right))
            basis(left) === basis(right) && left_domain == right_domain ||
                throw(ArgumentError("retained comparison basis or domain differs"))
            size(observe(row.error,Engine.absolute_error)) == size(observe(row.error,Engine.relative_error)) ==
                (length(reference.metadata.port_order),length(reference.metadata.port_order)) ||
                throw(DimensionMismatch("RMS matrices differ from output terminal identities"))
        end
        for key in definition.requested
            key in (:requests,:bands,:normalizations) && continue
            isequal(getproperty(definition.settings,key),getproperty(settings,key)) ||
                throw(ArgumentError("requested $key differs from retained analysis; request explicit reanalysis"))
        end
        selections=Tuple(key for key in definition.requested if key in (:requests,:bands,:normalizations))
        settings=merge(settings,NamedTuple{selections}(Tuple(getproperty(definition.settings,key) for key in selections)))
        pairs=unique((row.reference_index,row.candidate_index) for row in comparisons)
        if !isempty(selections)
            for request in settings.requests, band in settings.bands, normalization in settings.normalizations, pair in pairs
                any(row -> (row.reference_index,row.candidate_index)==pair && row.request==request &&
                    details(row.error).band==band && details(row.error).normalization==normalization,comparisons) ||
                    throw(ArgumentError("requested comparison was not retained; select an analysis or request explicit reanalysis"))
            end
            comparisons=filter(row -> row.request in settings.requests && details(row.error).band in settings.bands &&
                details(row.error).normalization in settings.normalizations,comparisons)
        end
    else
        comparisons=reference.result isa ObservationPublication || candidate.result isa ObservationPublication ?
            select(definition,reference,candidate) : select(definition,reference.result,candidate.result)
    end
    return (;reference,candidate,context,settings,comparisons,measurements)
end

"""Compare requests through their scientific owners, retaining point identity."""
function select(definition::BenchmarkTableDefinition,
        reference::Union{AbstractCoreResult,ParametricResult,AbstractUncertaintyResult,NamedTuple{(:result,:metadata)}},
        candidate::Union{AbstractCoreResult,ParametricResult,AbstractUncertaintyResult,NamedTuple{(:result,:metadata)}})
    settings=definition.settings
    comparisons=NamedTuple[]
    if settings.pairing !== nothing && reference isa Union{AbstractCoreResult,NamedTuple{(:result,:metadata)}} &&
            candidate isa Union{AbstractCoreResult,NamedTuple{(:result,:metadata)}}
        settings.pairing==((1,1),) || settings.pairing==[(1,1)] ||
            throw(ArgumentError("a scalar reference/candidate comparison permits only pairing=(1,1)"))
    end
    # Validate every requested product before calculating any comparison.
    for request in settings.requests
        observation_request(reference isa ParametricResult ? first(reference) : reference isa NamedTuple ? reference.result : reference,request)
        observation_request(candidate isa ParametricResult ? first(candidate) : candidate isa NamedTuple ? candidate.result : candidate,request)
    end
    for band in settings.bands, request in settings.requests, normalization in settings.normalizations
        controls=(;band,normalization,atol=settings.atol,fundamental=settings.fundamental,
            harmonics=settings.harmonics,unsupported=settings.unsupported)
        space_reference=reference isa Union{ParametricResult,AbstractUncertaintyResult}
        result=space_reference ?
            Engine.compare(reference,candidate,request;pairing=settings.pairing,controls...) :
            Engine.compare(reference,candidate,request;controls...)
        errors=result isa ParametricResult ? collect(result) :
            result isa AbstractVector ? result : (result,)
        identity=request_identity(request)
        quantity=Symbol(Units.symbol(request_quantity(request)))
        statistic=identity isa Tuple && first(identity) === UQ.statistics ?
            (last(identity) isa Base.Fix2{typeof(Statistics.quantile)} ?
                Symbol("q",lpad(string(round(Int,100last(identity).x)),2,'0')) :
                last(identity) === minimum ? :min : last(identity) === maximum ? :max : nameof(last(identity))) : :value
        for (index,error) in enumerate(errors)
            reference_index=settings.pairing === nothing ? 1 :
                first(only(filter(pair -> last(pair)==index,settings.pairing)))
            push!(comparisons,(request,quantity,statistic,reference_index,candidate_index=index,error))
        end
    end
    return comparisons
end
function validate(::BenchmarkTableDefinition,result::ObservationPublication,metadata::NamedTuple)
    basis(result)==metadata.basis || throw(ArgumentError("publication and operand basis differ"))
    n=length(metadata.port_order)
    all(payload -> size(payload.values)==(n,n,length(metadata.frequencies)),result) ||
        throw(DimensionMismatch("retained publication products do not match operand coordinates"))
    if haskey(result.columns,:frequency)
        contract=get(result.metadata.observation_columns,:frequency,nothing)
        scale=contract===nothing ? 1 : Units.scale_factor(contract.unit,Units.units(:base,:hertz))
        unique(result.columns.frequency).*scale==metadata.frequencies ||
            throw(ArgumentError("publication and operand frequency coordinates differ"))
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Expose comparison records, every matrix term and a labeled summary of per-term
maxima. Absolute RMS uses native units; relative RMS is displayed as percent.
Unavailable relative values retain their reasons and measured absolute values.
Term/maxima tables use scalar coordinates. Their unformatted maxima records
remain in `metadata(table.maxima, "comparison_records")` for the existing writer;
changing display columns therefore does not change the saved summary schema.
"""
function tabulate(definition::BenchmarkTableDefinition, source,
        published::NamedTuple{(:reference, :candidate, :context, :settings, :comparisons, :measurements)})
    calculations=DataFrame([(;role,case_id=string(published.context.case_id),
        path=get(operand.metadata,:path,missing),
        points=operand.result isa Union{ParametricResult,AbstractUncertaintyResult} ? length(operand.result) : 1,
        formulations=length(operand.metadata.formulation_sources),
        terminals=length(operand.metadata.port_order)) for (role,operand) in
        pairs((reference=published.reference,candidate=published.candidate))])
    comparisons=NamedTuple[]
    terms=NamedTuple[]
    summaries=NamedTuple[]
    formulations=NamedTuple[]
    selections=(published.reference.metadata.formulation_sources,published.candidate.metadata.formulation_sources)
    sources=Any[selections[1]...;selections[2]...]
    roles=vcat(fill(:reference,length(selections[1])),fill(:candidate,length(selections[2])))
    numbers=vcat(zeros(Int,length(selections[1])),collect(eachindex(selections[2])))
    labels=description(sources;roles,indices=numbers)
    formula_details=NamedTuple[]
    for (index,(role,source)) in enumerate(zip(roles,sources))
        number=role===:reference ? missing : numbers[index]
        push!(formulations,(role,formulation_index=number,label=labels[index]))
        ismissing(source) && continue
        for (scope,selected) in pairs((source isa Pair ? Tuple(source) : (source,))...)
            push!(formula_details,(role,formulation_index=number,label=labels[index],
                selection=description(scope,selected;compact=false),
                identifier=ismissing(formula_id(selected)) ? missing : string(formula_id(selected))))
        end
    end
    for (analysis,row) in enumerate(published.comparisons)
        error=row.error
        detail=details(error)
        ports=published.candidate.metadata.port_order
        left=select(published.reference.result,row.reference_index)
        right=select(published.candidate.result,row.candidate_index)
        candidate_axes=published.candidate.result isa ParametricResult ? published.candidate.result.axes : nothing
        nproblems=candidate_axes === nothing ?
            (published.candidate.result isa AbstractUncertaintyResult ? length(published.candidate.result) : 1) :
            length(candidate_axes.problems)
        problem_index=mod1(row.candidate_index,nproblems)
        formulation_index=cld(row.candidate_index,nproblems)
        unit=Units.native_unit(request_quantity(row.request),basis(right))
        absolute=detach(observe(error,Engine.absolute_error),1,definition.clip)
        relative=detach(observe(error,Engine.relative_error),100,definition.clip)
        statuses=get(detail,:status,fill(:not_recorded,size(absolute)))
        local_reasons=get(detail,:normalization_reason,fill(nothing,size(absolute)))
        reasons=map(reason -> reason === nothing ? get(detail,:reason,nothing) : reason,local_reasons)
        identity=(analysis, snapshot=get(row,:snapshot,"live"), benchmark=string(published.context.id), case_id=string(published.context.case_id),
            collection=string(published.context.collection), request=row.request, quantity=row.quantity, statistic=row.statistic,
            reference_point=row.reference_index, candidate_point=row.candidate_index,
            problem_index, formulation_index, band=detail.band, normalization=detail.normalization,
            absolute_unit=Units.label(unit),samples=detail.sample_count,
            resolution_revision=get(get(detail,:resolution,(;)),:revision,0),
            resolution_kind=get(get(detail,:resolution,(;)),:kind,:historical_unversioned),
            requested_bounds_Hz=get(detail,:requested_bounds,missing), actual_bounds_Hz=detail.actual_bounds)
        push!(comparisons,merge(identity,(absolute_rms=absolute,relative_rms_percent=relative,
            status=copy(statuses),reason=reasons,port_order=copy(ports),
            sample_indices=copy(get(detail,:indices,Int[])),tolerance=get(detail,:atol,missing),
            candidate_tolerance=get(detail,:candidate_atol,get(detail,:atol,missing)))))
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
            compared=length(ri),
            reasons=unique(filter(!isnothing,vec(reasons))))))
    end
    features=NamedTuple[]
    partitions=unique((;row.snapshot,row.benchmark,row.case_id,row.collection,row.problem_index,
        row.reference_point,row.request,row.quantity,row.statistic,row.normalization,row.absolute_unit) for row in summaries)
    for partition in partitions
        metrics=filter(row -> all(key -> getproperty(row,key)==getproperty(partition,key),keys(partition)),summaries)
        indices=unique(row.formulation_index for row in metrics)
        quantity_labels=description(sources;roles,indices=numbers,quantity=request_quantity(partition.request))
        selected_labels=quantity_labels[length(selections[1]).+indices]
        relative=DataFrame(formulation_index=indices,formula=selected_labels)
        absolute=copy(relative)
        for band in published.settings.bands
            name=band isa Symbol ? band : Symbol(string(band))
            for (frame,field) in ((relative,:maximum_relative_rms_percent),(absolute,:maximum_absolute_rms))
                frame[!,name]=map(indices) do index
                    matching=filter(row -> row.formulation_index==index && row.band==band,metrics)
                    isempty(matching) ? missing : getproperty(only(matching),field)
                end
            end
        end
        push!(features,merge(partition,(;relative,absolute)))
    end
    terms_table=DataFrame(terms)
    maxima_table=DataFrame(summaries)
    summary_table=tabulate(definition,source,summaries)
    # Keep structured comparison products in the audit table/publication. The
    # displayed term tables retain separately filterable physical coordinates.
    for frame in (terms_table,maxima_table,summary_table)
        frame[!,:band]=[value isa Symbol ? value : string(value) for value in frame.band]
        for (field,prefix) in ((:requested_bounds_Hz,"requested"),(:actual_bounds_Hz,"actual"))
            bounds=frame[!,field]
            frame[!,Symbol(prefix*"_lower_Hz")]=[value===nothing || ismissing(value) ? missing : first(value) for value in bounds]
            frame[!,Symbol(prefix*"_upper_Hz")]=[value===nothing || ismissing(value) ? missing : last(value) for value in bounds]
        end
        select!(frame,Not([:requested_bounds_Hz,:actual_bounds_Hz]))
    end
    for frame in (terms_table,maxima_table)
        names_by_quantity=Dict(quantity => description(sources;roles,indices=numbers,quantity)
            for quantity in unique(request_quantity.(frame.request)))
        frame[!,:method]=[names_by_quantity[request_quantity(row.request)][length(selections[1])+row.formulation_index] for row in eachrow(frame)]
        reference_count=published.reference.result isa ParametricResult ? length(published.reference.result.axes.problems) :
            published.reference.result isa AbstractUncertaintyResult ? length(published.reference.result) : 1
        frame[!,:reference_method]=[names_by_quantity[request_quantity(row.request)][cld(row.reference_point,reference_count)] for row in eachrow(frame)]
        select!(frame,Not(:request))
    end
    terms_table[!,:reason]=[value===nothing ? missing : value for value in terms_table.reason]
    for field in (:absolute_term,:relative_term)
        maxima_table[!,Symbol(string(field)*"_response")]=[ismissing(value) ? missing : first(value) for value in maxima_table[!,field]]
        maxima_table[!,Symbol(string(field)*"_excitation")]=[ismissing(value) ? missing : last(value) for value in maxima_table[!,field]]
    end
    maxima_table[!,:reasons]=[join(string.(reasons),"; ") for reasons in maxima_table.reasons]
    select!(maxima_table,Not([:absolute_term,:relative_term]))
    metadata!(maxima_table,"comparison_records",summaries;style=:note)
    method_labels=(reference=join(labels[1:length(selections[1])]," / "),
        candidate=join(labels[length(selections[1])+1:end]," / "))
    return merge((calculations,formulations=DataFrame(formulations),formula_details=DataFrame(formula_details),comparisons=DataFrame(comparisons),
        terms=terms_table,maxima=maxima_table,summary=summary_table,features),
        tabulate(definition,published.measurements;labels=method_labels,
            indices=(reference=missing,candidate=length(selections[2])==1 ? 1 : missing)),
        tabulate(definition,(reference=published.reference,candidate=published.candidate);labels=method_labels))
end

function illustrate(definition::BenchmarkTableDefinition, source,
        published::NamedTuple{(:reference, :candidate, :context, :settings, :comparisons, :measurements)}, table)
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
function validate(definition::BenchmarkTableDefinition,result::AbstractUncertaintyResult,metadata::NamedTuple)
    isempty(result) && throw(ArgumentError("uncertainty result cannot be empty"))
    foreach(value -> validate(definition,value,metadata),result)
    return nothing
end
