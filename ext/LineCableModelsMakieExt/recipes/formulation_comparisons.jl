import LineCableModels.ReportBuilder: ReportArtifact, select

# Identify requested defaults structurally. Fixed choices shared by the whole
# catalogue do not prevent its varying formula choices from having a default.
function _addon_default_formulations(records)
    identifiers = map(records) do record
        fields = Dict{Tuple,Any}()
        function visit(value, path)
            value isa NamedTuple || return
            if haskey(value, :identifier)
                fields[path] = value.identifier
            else
                for (key, child) in pairs(value)
                    visit(child, (path..., key))
                end
            end
        end
        visit(get(record, :requested, (;)), ())
        fields
    end
    paths = unique([path for fields in identifiers for path in keys(fields)])
    varying = filter(path -> !all(fields -> get(fields, path, nothing) ==
        get(first(identifiers), path, nothing), identifiers), paths)
    considered = isempty(varying) ? paths : varying
    return [!isempty(considered) && all(path -> get(fields, path, nothing) === :default,
        considered) for fields in identifiers]
end

function _addon_candidate_style_indices(defaults)
    first_default = findfirst(defaults)
    return [2 * (first_default === nothing ? index : index == first_default ? 1 :
        index == 1 ? first_default : index) for index in eachindex(defaults)]
end

# Project only the requested equation choices, never the curves or their stable
# catalogue indices. Shared physical/numerical controls remain in the labels.
function _addon_formulation_labels(records, family; kwargs...)
    omitted = family === Val(:series) ?
        (:earth_admittance, :insulation_admittance, :semicon_admittance) :
        (:earth_impedance, :internal_impedance, :insulation_impedance, :pipe_impedance)
    projected = map(records) do record
        haskey(record, :requested) || return record
        requested = (; (key => value for (key, value) in pairs(record.requested)
            if key ∉ omitted)...)
        merge(record, (; requested))
    end
    return LineCableModels.description(projected; kwargs...)
end

"""
Plot completed formulation results for an explicitly selected problem. All
formulations are overlaid by default; filtering retains original labels/colors.
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
    records=[value isa NamedTuple ? value : NamedTuple(value) for value in results.axes.formulations]
    labels=LineCableModels.description(records)
    family_labels=Dict(family => _addon_formulation_labels(records,family)
        for family in (Val(:series),Val(:shunt)))
    defaults=_addon_default_formulations(records)
    catalogue_styles=_addon_candidate_style_indices(defaults)
    built=LineCableModels.UIPlot[]
    for p in problems
        sources=Tuple(results[p,index] for index in indices)
        names=Tuple(labels[index] for index in indices)
        page_labels=Dict(family => Tuple(values[index] for index in indices)
            for (family,values) in family_labels)
        styles=Tuple(catalogue_styles[index] for index in indices)
        roles=Tuple(defaults[index] ? :default : :alternative for index in indices)
        if reference !== nothing
            reference isa LineCableModels.LineParameters || throw(ArgumentError("reference must be a scalar LineParameters result"))
            sources=(reference,sources...)
            reference_record=get(LineCableModels.details(reference),:formulations,(;))
            names=(only(LineCableModels.description([reference_record];prefix="Reference F")),names...)
            page_labels=Dict(family => (only(_addon_formulation_labels([reference_record],family;
                prefix="Reference F")),values...) for (family,values) in page_labels)
            styles=(1,styles...)
            roles=(:reference,roles...)
        end
        all(value -> value isa LineCableModels.LineParameters,sources) || throw(ArgumentError("matrix-curve overlays require line-parameter results"))
        normalized=_line_plot_ydata(first(sources),selected_ydata)
        pages=_addon_line_pages(sources; ydata=normalized, series_labels=series_labels === nothing ? names : series_labels,
            series_family_labels=series_labels === nothing ? page_labels : nothing,
            series_indices=styles,xscale=_scale_symbol(xscale),yscale=_scale_symbol(yscale),clip,
            series_defaults=_addon_comparison_styles(styles,roles,2length(records)),
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
function plot(published::NamedTuple{(:reference,:candidate,:context,:settings,:comparisons)},
        selection=nothing; ydata=nothing, problem=nothing, formulations=nothing,
        pair=nothing, band=nothing, series_labels=nothing, xscale=:log10, yscale=:linear,
        clip=true,atol=published.settings.atol,
        legend_position=:bottom,legend_overflow=:show_all,legend_attributes=(;),kwargs...)
    selected_ydata=_plot_ydata(selection,ydata,
        (LineCableModels.Z,LineCableModels.Y))
    current_resolution = all(row -> get(get(row.error.details,:resolution,(;)),:revision,0) ==
        LineCableModels.Engine.OBSERVABLE_RESOLUTION_REVISION, published.comparisons)
    current_resolution || @warn "Historical comparison semantics: curves use current observation resolution; retained RMS is unchanged. Request explicit reanalysis for a current report."
    isequal(atol,published.settings.atol) || @warn "Plot resolution override differs from retained comparison controls; retained RMS is unchanged." atol
    candidate=published.candidate.result
    reference=published.reference.result
    isspace=candidate isa LineCableModels.ParametricResult
    nproblems=isspace ? length(candidate.axes.problems) : 1
    problem === nothing && nproblems != 1 && throw(ArgumentError("select problem explicitly before plotting multiple problems"))
    problems=problem === nothing ? [1] : problem isa Integer ? [problem] : collect(problem)
    all(index -> index isa Integer && 1 <= index <= nproblems,problems) || throw(ArgumentError("invalid problem selection"))
    count=isspace ? length(candidate.axes.formulations) : 1
    indices=formulations === nothing ? collect(1:count) : formulations isa Integer ? [formulations] : collect(formulations)
    !isempty(indices) && allunique(indices) && all(index -> index isa Integer && 1 <= index <= count,indices) ||
        throw(ArgumentError("invalid formulation selection"))
    declared=unique([(row.reference_index,row.candidate_index) for row in published.comparisons])
    pair === nothing || pair in declared || throw(ArgumentError("pair must identify an explicitly saved reference/candidate comparison"))
    raw=isspace ? candidate.axes.formulations : [published.candidate.metadata.formulation]
    records=[value isa NamedTuple ? value : NamedTuple(value) for value in raw]
    reference_records=reference isa LineCableModels.ParametricResult ?
        [value isa NamedTuple ? value : NamedTuple(value) for value in reference.axes.formulations] :
        [published.reference.metadata.formulation]
    reference_problems=reference isa LineCableModels.ParametricResult ? length(reference.axes.problems) : 1
    reference_labels=LineCableModels.description(reference_records;prefix="Reference F")
    labels=LineCableModels.description(records)
    family_labels=Dict(family => (
        reference=_addon_formulation_labels(reference_records,family;prefix="Reference F"),
        candidate=_addon_formulation_labels(records,family))
        for family in (Val(:series),Val(:shunt)))
    defaults=_addon_default_formulations(records)
    catalogue_styles=_addon_candidate_style_indices(defaults)
    built=LineCableModels.UIPlot[]
    for p in problems
        points=[p+(f-1)*nproblems for f in indices]
        selected=[entry for entry in declared if last(entry) in points && (pair === nothing || entry==pair)]
        isempty(selected) && throw(ArgumentError("no retained comparisons match the selection"))
        refs=unique(first.(selected))
        candidates=unique(last.(selected))
        # Equal numerical curves are distinct declared selections, never set elements.
        sources=Tuple(vcat([select(reference,i) for i in refs],[select(candidate,i) for i in candidates]))
        all(value -> value isa LineCableModels.LineParameters,sources) || throw(ArgumentError("matrix-curve overlays require line-parameter operands; moment errors remain available through report"))
        names=Tuple(vcat([reference_labels[cld(i,reference_problems)] for i in refs],[labels[cld(i,nproblems)] for i in candidates]))
        page_labels=Dict(family => Tuple(vcat(
            [values.reference[cld(i,reference_problems)] for i in refs],
            [values.candidate[cld(i,nproblems)] for i in candidates]))
            for (family,values) in family_labels)
        styles=Tuple(vcat([2cld(i,reference_problems)-1 for i in refs],
            [catalogue_styles[cld(i,nproblems)] for i in candidates]))
        roles=Tuple(vcat(fill(:reference,length(refs)),
            [defaults[cld(i,nproblems)] ? :default : :alternative for i in candidates]))
        if band !== nothing
            sources=map(sources, vcat([(:reference,i) for i in refs],[(:candidate,i) for i in candidates])) do value,entry
                role,index=entry
                rows=filter(row -> row.error.details.band == band &&
                    (role === :reference ? row.reference_index : row.candidate_index)==index,published.comparisons)
                isempty(rows) && throw(ArgumentError("band was not retained; request explicit reanalysis before plotting it"))
                samples=first(rows).error.details.indices
                isempty(samples) && throw(ArgumentError("the selected band has no retained samples"))
                all(row -> row.error.details.indices == samples,rows) || throw(ArgumentError("selected band has conflicting retained sample coordinates"))
                value[samples]
            end |> Tuple
        end
        normalized=_line_plot_ydata(first(sources),selected_ydata)
        pages=_addon_line_pages(sources;ydata=normalized,series_labels=series_labels === nothing ? names : series_labels,
            series_family_labels=series_labels === nothing ? page_labels : nothing,
            series_indices=styles,xscale=_scale_symbol(xscale),yscale=_scale_symbol(yscale),clip,atol,
            series_defaults=_addon_comparison_styles(styles,roles,2max(length(records),length(reference_records))),
            legend_position,legend_overflow,legend_attributes,
            signed_ylog=true,kwargs...)
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
