import LineCableModels.ReportBuilder: ReportArtifact, select

# Styles consume owner-scoped identifiers, not serialized field layouts.
function _addon_default_formulations(sources;quantity=nothing)
    return map(sources) do source
        ismissing(source) && return false
        identifiers=[LineCableModels.formula_id(value) for (scope,value) in
            pairs((source isa Pair ? Tuple(source) : (source,))...;quantity) if !isempty(last(scope))]
        !isempty(identifiers) && all(id -> id===:default,identifiers)
    end
end

"""
Plot completed formulation results for an explicitly selected problem. Unique
quantity-relevant formulations are overlaid; filtering retains source order/colors.
A scalar reference is drawn once. No solve or comparison is performed.
"""
function plot(results::LineCableModels.ParametricResult, selection=nothing;
        ydata=nothing,
        problem=nothing, formulations=nothing, reference=nothing,
        series_labels=nothing, xscale=:log10, yscale=:linear, clip=true,
        legend_position=:bottom, legend_overflow=:show_all, legend_attributes=(;), kwargs...)
    selected_ydata=_plot_ydata(selection,ydata,
        (LineCableModels.Z,LineCableModels.Y))
    isempty(results.axes) && throw(ArgumentError("plotting requires retained problem/formulation axes"))
    count=length(results.axes.problems)
    problem === nothing && count != 1 && throw(ArgumentError("select problem explicitly before plotting multiple problems"))
    problems=problem === nothing ? [1] : problem isa Integer ? [problem] : collect(problem)
    all(index -> index isa Integer && 1 <= index <= count,problems) || throw(ArgumentError("invalid problem selection"))
    indices=formulations === nothing ? collect(eachindex(results.axes.formulations)) :
        formulations isa Integer ? [formulations] : collect(formulations)
    !isempty(indices) && allunique(indices) && all(index -> index isa Integer &&
        1 <= index <= length(results.axes.formulations),indices) || throw(ArgumentError("invalid formulation selection"))
    records=[ImportExport.deserialize_value(Val(:formulation),value) for value in results.axes.formulations]
    labels=LineCableModels.description(records)
    family_labels=Dict(family => LineCableModels.description(records;quantity)
        for (family,quantity) in ((Val(:series),LineCableModels.Z),(Val(:shunt),LineCableModels.Y)))
    built=LineCableModels.UIPlot[]
    for p in problems
        sources=Tuple(results[p,index] for index in indices)
        names=Tuple(labels[index] for index in indices)
        page_labels=Dict(family => Tuple(values[index] for index in indices)
            for (family,values) in family_labels)
        styles=Tuple(index+1 for index in indices)
        displayed_records=records[indices]
        roles=fill(:candidate,length(indices))
        if reference !== nothing
            reference isa LineCableModels.LineParameters || throw(ArgumentError("reference must be a scalar LineParameters result"))
            sources=(reference,sources...)
            reference_record=ImportExport.deserialize_value(Val(:formulation),
                get(LineCableModels.details(reference),:formulations,(;)))
            all_sources=Any[reference_record;records]
            label_roles=[:reference;fill(:candidate,length(records))]
            all_labels=LineCableModels.description(all_sources;roles=label_roles)
            names=Tuple(all_labels[[1;indices.+1]])
            page_labels=Dict(family => Tuple(LineCableModels.description(all_sources;
                roles=label_roles,quantity)[[1;indices.+1]])
                for (family,quantity) in ((Val(:series),LineCableModels.Z),(Val(:shunt),LineCableModels.Y)))
            styles=(1,styles...)
            roles=[:reference;roles]
            displayed_records=Any[reference_record;displayed_records]
        end
        all(value -> value isa LineCableModels.LineParameters,sources) || throw(ArgumentError("matrix-curve overlays require line-parameter results"))
        normalized=_line_plot_ydata(first(sources),selected_ydata)
        pages=_addon_line_pages(sources; ydata=normalized, series_labels=series_labels === nothing ? names : series_labels,
            series_family_labels=series_labels === nothing ? page_labels : nothing,
            series_indices=styles,xscale=xscale,yscale=yscale,clip,
            formulation_sources=displayed_records,formulation_roles=roles,
            legend_position,legend_overflow,legend_attributes,
            signed_ylog=true,kwargs...)
        for page in (pages isa LineCableModels.UIPlot ? (pages,) : pages)
            page.addon_state=merge(page.addon_state,(formulations=(problem=p,indices=copy(indices),records=records[indices],reference=reference === nothing ? nothing : LineCableModels.details(reference)),))
            push!(built,page)
        end
    end
    return length(built)==1 ? only(built) : built
end

"""Plot selected completed benchmark data; report construction is never repeated."""
function plot(artifact::ReportArtifact, selection=nothing; ydata=nothing, kwargs...)
    selected_ydata=_plot_ydata(selection,ydata,
        (LineCableModels.Z,LineCableModels.Y))
    return plot(artifact.published,selected_ydata;kwargs...)
end

"""Overlay retained reference/candidate comparisons for explicitly selected problems."""
function plot(published::NamedTuple{(:reference,:candidate,:context,:settings,:comparisons,:measurements)},
        selection=nothing; ydata=nothing, problem=nothing, formulations=nothing,
        pair=nothing, band=nothing, series_labels=nothing, xscale=:log10, yscale=:linear,
        clip=true,atol=published.settings.atol,
        legend_position=:bottom,legend_overflow=:show_all,legend_attributes=(;),kwargs...)
    selected_ydata=_plot_ydata(selection,ydata,
        (LineCableModels.Z,LineCableModels.Y))
    current_resolution = all(row -> get(get(LineCableModels.details(row.error),:resolution,(;)),:revision,0) ==
        Engine.OBSERVABLE_RESOLUTION_REVISION, published.comparisons)
    current_resolution || @warn "Historical comparison semantics: curves use current observation resolution; retained RMS is unchanged. Request explicit reanalysis for a current report."
    isequal(atol,published.settings.atol) || @warn "Plot resolution override differs from retained comparison controls; retained RMS is unchanged." atol
    candidate=published.candidate.result
    reference=published.reference.result
    isspace=candidate isa LineCableModels.ParametricResult
    isuQ=candidate isa Union{LineCableModels.AbstractUncertaintyResult,ObservationPublication}
    desired = selected_ydata isa Function ||
        (selected_ydata isa Tuple && !isempty(selected_ydata) && first(selected_ydata) isa Function &&
            (first(selected_ydata) === LineCableModels.statistics || any(item -> !(item isa Function) || item isa Colon,selected_ydata))) ?
        (selected_ydata,) : selected_ydata
    statistical = map(desired) do item
        identity = request_identity(item)
        identity isa Tuple && !isempty(identity) && first(identity) === LineCableModels.statistics
    end
    any(statistical) && !all(statistical) && throw(ArgumentError(
        "plot ordinary quantities and explicit statistics in separate calls"))
    statistical_view = isuQ && !isempty(statistical) && all(statistical)
    isuQ && !statistical_view &&
        (reference isa ObservationPublication || candidate isa ObservationPublication) &&
        throw(ArgumentError("mean ± std overlays require uncertainty-bearing core results; select explicit retained statistics for this publication"))
    nproblems=isspace ? length(candidate.axes.problems) : candidate isa LineCableModels.AbstractUncertaintyResult ? length(candidate) : 1
    problem === nothing && nproblems != 1 && throw(ArgumentError("select problem explicitly before plotting multiple problems"))
    problems=problem === nothing ? [1] : problem isa Integer ? [problem] : collect(problem)
    all(index -> index isa Integer && 1 <= index <= nproblems,problems) || throw(ArgumentError("invalid problem selection"))
    count=isspace ? length(candidate.axes.formulations) : 1
    indices=formulations === nothing ? collect(1:count) : formulations isa Integer ? [formulations] : collect(formulations)
    !isempty(indices) && allunique(indices) && all(index -> index isa Integer && 1 <= index <= count,indices) ||
        throw(ArgumentError("invalid formulation selection"))
    declared=unique([(row.reference_index,row.candidate_index) for row in published.comparisons])
    pair === nothing || pair in declared || throw(ArgumentError("pair must identify an explicitly saved reference/candidate comparison"))
    records=published.candidate.metadata.formulation_sources
    reference_records=published.reference.metadata.formulation_sources
    reference_problems=reference isa LineCableModels.ParametricResult ? length(reference.axes.problems) :
        reference isa LineCableModels.AbstractUncertaintyResult ? length(reference) : 1
    all_sources=Any[reference_records...;records...]
    label_roles=vcat(fill(:reference,length(reference_records)),fill(:candidate,length(records)))
    combined=LineCableModels.description(all_sources;roles=label_roles)
    reference_labels=combined[1:length(reference_records)]
    labels=combined[length(reference_records)+1:end]
    family_labels=Dict(family => let
        values=LineCableModels.description(all_sources;roles=label_roles,quantity)
        (reference=values[1:length(reference_records)],candidate=values[length(reference_records)+1:end])
    end for (family,quantity) in ((Val(:series),LineCableModels.Z),(Val(:shunt),LineCableModels.Y)))
    built=LineCableModels.UIPlot[]
    for p in problems
        points=[p+(f-1)*nproblems for f in indices]
        selected=[entry for entry in declared if last(entry) in points && (pair === nothing || entry==pair)]
        isempty(selected) && throw(ArgumentError("no retained comparisons match the selection"))
        refs=unique(first.(selected))
        candidates=[point for point in points if any(entry -> last(entry)==point,selected)]
        # Raw points remain intact; each quantity page selects its unique formulas.
        read_core = isuQ && !statistical_view ? LineCableModels.uncertain : select
        sources=Tuple(vcat([read_core(reference,i) for i in refs],[read_core(candidate,i) for i in candidates]))
        all(value -> value isa Union{LineCableModels.LineParameters,ObservationPublication},sources) || throw(ArgumentError("matrix-curve overlays require retained matrix coordinates"))
        names=Tuple(vcat([reference_labels[cld(i,reference_problems)] for i in refs],[labels[cld(i,nproblems)] for i in candidates]))
        page_labels=Dict(family => Tuple(vcat(
            [values.reference[cld(i,reference_problems)] for i in refs],
            [values.candidate[cld(i,nproblems)] for i in candidates]))
            for (family,values) in family_labels)
        styles=Tuple(vcat([cld(i,reference_problems) for i in refs],
            [length(reference_records)+cld(i,nproblems) for i in candidates]))
        roles=Tuple(vcat(fill(:reference,length(refs)),
            fill(:candidate,length(candidates))))
        displayed_records=Any[[reference_records[cld(i,reference_problems)] for i in refs]...;
            [records[cld(i,nproblems)] for i in candidates]...]
        if band !== nothing && !statistical_view
            sources=map(sources, vcat([(:reference,i) for i in refs],[(:candidate,i) for i in candidates])) do value,entry
                role,index=entry
                rows=filter(row -> LineCableModels.details(row.error).band == band &&
                    (role === :reference ? row.reference_index : row.candidate_index)==index,published.comparisons)
                isempty(rows) && throw(ArgumentError("band was not retained; request explicit reanalysis before plotting it"))
                samples=LineCableModels.details(first(rows).error).indices
                isempty(samples) && throw(ArgumentError("the selected band has no retained samples"))
                all(row -> LineCableModels.details(row.error).indices == samples,rows) || throw(ArgumentError("selected band has conflicting retained sample coordinates"))
                value[samples]
            end |> Tuple
        end
        prepared=nothing
        normalized=if statistical_view
            reference isa Union{LineCableModels.AbstractUncertaintyResult,ObservationPublication} || throw(ArgumentError("UQ overlays require explicit compatible statistical operands"))
            requests=Tuple(request_identity(item)==request ? item : request
                for item in desired for request in published.settings.requests if
                request_identity(item)==request)
            isempty(requests) && throw(ArgumentError("no retained statistical requests match ydata"))
            allunique(requests) || throw(ArgumentError("duplicate statistical plot requests"))
            rows=filter(row -> LineCableModels.details(row.error).band==band,published.comparisons)
            band===nothing || !isempty(rows) || throw(ArgumentError("band was not retained; request explicit reanalysis before plotting it"))
            samples=band===nothing ? nothing : LineCableModels.details(first(rows).error).indices
            band===nothing || all(row -> LineCableModels.details(row.error).indices==samples,rows) ||
                throw(ArgumentError("selected band has conflicting retained sample coordinates"))
            options=(; (key=>value for (key,value) in kwargs if key in (:freq_unit,:length_unit,:quantity_units))...)
            prepared=Tuple(vcat([_prepare_line_observations(reference;point=i,ydata=requests,sample_indices=samples,clip,atol,options...) for i in refs],
                [_prepare_line_observations(candidate;point=i,ydata=requests,sample_indices=samples,clip,atol,options...) for i in candidates]))
            requests
        else
            _line_plot_ydata(first(sources),selected_ydata)
        end
        pages=_addon_line_pages(sources;publications=prepared,ydata=normalized,series_labels=series_labels === nothing ? names : series_labels,
            series_family_labels=series_labels === nothing ? page_labels : nothing,
            series_indices=styles,xscale=xscale,yscale=yscale,clip,atol,
            formulation_sources=displayed_records,formulation_roles=roles,
            legend_position,legend_overflow,legend_attributes,
            title_prefix=get(published.context,:id,nothing),signed_ylog=true,kwargs...)
        for page in (pages isa LineCableModels.UIPlot ? (pages,) : pages)
            page.addon_state=merge(page.addon_state,(
                resolution=(atol,clip,current_comparison=current_resolution,
                    display_override=!isequal(atol,published.settings.atol)),
                formulations=(problem=p,indices=copy(indices),records=records[indices],references=reference_records),))
            push!(built,page)
        end
    end
    return length(built)==1 ? only(built) : built
end
