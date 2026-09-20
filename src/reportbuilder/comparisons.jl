function _selected_errors(definition,point,reference)
    rows=point.errors
    isempty(rows) && throw(ArgumentError("this observation has no completed comparisons"))
    if reference!==nothing
        rows=filter(row -> row.reference_id==reference.gridpoint.id,rows)
    end
    settings=definition.settings
    for key in definition.requested
        if key===:requests
            rows=filter(row -> row.request in settings.requests,rows)
        elseif key===:bands
            rows=filter(row -> row.band in settings.bands,rows)
        elseif key===:normalizations
            rows=filter(row -> row.normalization in settings.normalizations,rows)
        elseif key!==:pairing
            all(row -> isequal(get(row.settings,key===:atol ? :requested_atol : key,nothing),getproperty(settings,key)),rows) ||
                throw(ArgumentError("requested $key differs from the retained comparison; run compare explicitly"))
        end
    end
    isempty(rows) && throw(ArgumentError("requested comparison was not retained"))
    requests=:requests in definition.requested ? settings.requests : unique(row.request for row in rows)
    bands=:bands in definition.requested ? settings.bands : unique(row.band for row in rows)
    normalizations=:normalizations in definition.requested ? settings.normalizations : unique(row.normalization for row in rows)
    for request in requests,band in bands,normalization in normalizations
        any(row -> row.request==request && isequal(row.band,band) && row.normalization==normalization,rows) ||
            throw(ArgumentError("some requested comparison products were not retained"))
    end
    return rows
end

function _point_information(points,role)
    labels=Grammar.observation_labels(points)
    calculations=NamedTuple[]
    formulations=NamedTuple[]
    formula_details=NamedTuple[]
    for (index,point) in enumerate(points)
        id=point.gridpoint.id
        push!(calculations,(role,identity=id,label=labels[index],physical_inputs=point.gridpoint.inputs))
        push!(formulations,(role,identity=id,label=labels[index],selections=point.gridpoint.formulations))
        for (key,value) in pairs(something(point.gridpoint.formulations,(;)))
            push!(formula_details,(role,identity=id,selection=key,value))
        end
    end
    return (;calculations,formulations,formula_details)
end

"""
$(TYPEDSIGNATURES)

Tabulate completed comparisons, independent maxima, sampling information, and
recorded performance evidence. Every candidate remains represented. Display
features consume the same observation-side groups as plots; no source access,
comparison, sampling estimate, or timing measurement occurs here.
"""
function tabulate(definition::BenchmarkTableDefinition,
        observed::Union{ObservedResult,AbstractVector{<:ObservedResult}};reference=nothing)
    points=collect(_observed_points(observed))
    labels=Grammar.observation_labels(points)
    rows=NamedTuple[];terms=NamedTuple[];maxima=NamedTuple[]
    for (index,point) in enumerate(points),error in _selected_errors(definition,point,reference)
        id=point.gridpoint.id
        settings=error.settings
        identity=(candidate_id=error.candidate_id,reference_id=error.reference_id,
            problem_index=id.problem_index,formulation_index=id.formulation_index,
            candidate_point=index,method=labels[index],request=error.request,
            quantity=Symbol(Units.symbol(error.quantity)),statistic=error.statistic,
            band=error.band,normalization=error.normalization,
            absolute_unit=Units.label(error.absolute_unit),samples=settings.sample_count,
            requested_bounds_Hz=get(settings,:requested_bounds,nothing),actual_bounds_Hz=settings.actual_bounds)
        absolute=error.absolute
        relative=error.relative.*100
        status=get(settings,:status,fill(:not_recorded,size(absolute)))
        reasons=get(settings,:normalization_reason,fill(get(settings,:reason,nothing),size(absolute)))
        push!(rows,merge(identity,(absolute_rms=Grammar.detach(absolute),relative_rms_percent=relative,
            status=Grammar.detach(status),reason=Grammar.detach(reasons),settings)))
        for i in axes(absolute,1),j in axes(absolute,2)
            push!(terms,merge(identity,(row=i,column=j,response=error.coordinates[i],excitation=error.coordinates[j],
                absolute_rms=absolute[i,j],relative_rms_percent=relative[i,j],status=status[i,j],reason=reasons[i,j])))
        end
        push!(maxima,merge(identity,(maximum_absolute_rms=error.maxima.absolute.value,
            absolute_term=error.maxima.absolute.index,maximum_relative_rms_percent=error.maxima.relative.value*100,
            relative_term=error.maxima.relative.index,
            absolute_term_response=ismissing(error.maxima.absolute.index) ? missing : error.coordinates[error.maxima.absolute.index[1]],
            absolute_term_excitation=ismissing(error.maxima.absolute.index) ? missing : error.coordinates[error.maxima.absolute.index[2]],
            relative_term_response=ismissing(error.maxima.relative.index) ? missing : error.coordinates[error.maxima.relative.index[1]],
            relative_term_excitation=ismissing(error.maxima.relative.index) ? missing : error.coordinates[error.maxima.relative.index[2]],
            unavailable=count(ismissing,relative),
            term_count=length(relative),compared=count(!ismissing,relative),
            reasons=unique(filter(!isnothing,vec(reasons))))))
    end
    features=NamedTuple[]
    groups=NamedTuple[]
    partitions=unique((row.request,row.normalization,row.reference_id,row.problem_index) for row in rows)
    for (request,normalization,reference_id,problem_index) in partitions
        selected=filter(row -> row.request==request && row.normalization==normalization && row.reference_id==reference_id && row.problem_index==problem_index,maxima)
        bands=unique(row.band for row in selected)
        active=unique(row.candidate_point for row in selected)
        display_groups=[Grammar.observation_groups(points[active];request,band,normalization,reference=reference_id) for band in bands]
        representatives=sort(unique(vcat(([active[group.representative] for group in entries] for entries in display_groups)...)))
        push!(groups,(request,normalization,reference_id,bands,groups=display_groups,points=active))
        quantity_labels=Grammar.observation_labels(points;request)
        absolute=DataFrame(formula=quantity_labels[representatives]);relative=copy(absolute)
        for band in bands
            name=band isa Symbol ? band : Symbol(string(band))
            for (frame,field) in ((absolute,:maximum_absolute_rms),(relative,:maximum_relative_rms_percent))
                frame[!,name]=[begin
                    matches=filter(row -> row.candidate_point==index && isequal(row.band,band),selected)
                    isempty(matches) ? missing : getproperty(only(matches),field)
                end for index in representatives]
            end
        end
        push!(features,(request,quantity=first(selected).quantity,statistic=first(selected).statistic,
            normalization,reference_id,problem_index,absolute_unit=first(selected).absolute_unit,absolute,relative))
    end
    candidate_info=_point_information(points,:candidate)
    reference_info=reference===nothing ? (calculations=NamedTuple[],formulations=NamedTuple[],formula_details=NamedTuple[]) :
        _point_information([reference],:reference)
    info=map((a,b) -> DataFrame(vcat(a,b)),candidate_info,reference_info)
    timing=_timing_tables(points,reference)
    sampling=_sampling_tables(points,reference)
    coverage=DataFrame([(candidate_source=string(row.candidate_id.source_id),reference_source=string(row.reference_id.source_id),
        point=row.problem_index,formulation_index=row.formulation_index,band=row.band isa Tuple ? string(row.band) : row.band,
        frequency_count=row.samples,first_Hz=row.actual_bounds_Hz===nothing ? missing : Grammar.nominal(first(row.actual_bounds_Hz)),
        last_Hz=row.actual_bounds_Hz===nothing ? missing : Grammar.nominal(last(row.actual_bounds_Hz))) for row in maxima])
    coverage=unique(coverage)
    overview=(coverage,execution=timing.execution,performance=timing.performance,
        timing_ratio=timing.performance_comparison,source_timings=timing.source_timings,
        sampling=sampling.sampling,cdf_precision=isempty(sampling.sampling) ? DataFrame() : DataFrames.select(sampling.sampling,
            :role,:point,:method,:trials,:confidence,:marginal_count,:target_cdf,:cdf_bound,:target_supported,:scope))
    maxima_table=DataFrame(maxima)
    metadata!(maxima_table,"comparison_records",Grammar.detach(maxima);style=:note)
    return merge(info,(quantities=map(tabulate,points),comparisons=DataFrame(rows),terms=DataFrame(terms),
        maxima=maxima_table,summary=copy(maxima_table),features,groups),timing,sampling,(;overview))
end

function illustrate(definition::BenchmarkTableDefinition,observed,tables;reference=nothing)
    illustration=definition.illustration
    (illustration===nothing || illustration===false) && return nothing
    callable=illustration===true ? PlotBuilder.plot : illustration
    return callable(observed;reference,definition.plot_options...)
end

"""
$(TYPEDSIGNATURES)

Convenience for raw benchmark operands. Explicit comparison completes before
observation construction. Already completed products supplied as `comparisons`
are joined by original identities. The reference remains a separate observation.
"""
function report(definition::BenchmarkTableDefinition,source::NamedTuple;requests::Tuple=(),observation_options::NamedTuple=(;))
    all(key -> haskey(source,key),(:reference,:candidate)) || throw(ArgumentError("benchmark operands require reference and candidate"))
    reference=source.reference isa NamedTuple && haskey(source.reference,:result) ? source.reference.result : source.reference
    candidate=source.candidate isa NamedTuple && haskey(source.candidate,:result) ? source.candidate.result : source.candidate
    candidate isa Union{ObservedResult,AbstractVector{<:ObservedResult}} && return report(definition,candidate;reference)
    settings=definition.settings
    completed=haskey(source,:comparisons) ? source.comparisons : Engine.compare(reference,candidate,collect(settings.requests);
        bands=settings.bands,normalizations=settings.normalizations,pairing=settings.pairing,
        atol=settings.atol,fundamental=settings.fundamental,harmonics=settings.harmonics,unsupported=settings.unsupported)
    timings=(measurements=get(source,:measurements,nothing),context=get(source,:context,(;)))
    candidate isa AbstractUncertaintyResult && isempty(requests) && (requests=settings.requests)
    observed=observables(candidate,requests;comparisons=completed,timings,clip=definition.clip,atol=settings.atol,complete_pairs=true,observation_options...)
    observed_reference=reference isa AbstractUncertaintyResult ?
        (length(reference)==1 ? ObservedResult(reference,1,requests;clip=definition.clip,atol=settings.atol,complete_pairs=true,observation_options...) :
            throw(ArgumentError("a benchmark report retains one separate reference point"))) :
        ObservedResult(reference,requests;clip=definition.clip,atol=settings.atol,complete_pairs=true,observation_options...)
    return report(definition,observed;reference=observed_reference)
end

description(::BenchmarkTableDefinition,band) = band isa Symbol ? string(band) : string(first(band),"–",last(band)," Hz")
