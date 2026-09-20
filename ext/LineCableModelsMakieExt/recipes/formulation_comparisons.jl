import LineCableModels.ReportBuilder: ReportArtifact

"""Construct observations, then overlay selected original study points."""
function plot(results::LineCableModels.ParametricResult,selection=nothing;
        ydata=nothing,problem=nothing,formulations=nothing,reference=nothing,
        freq_unit=:base,length_unit=:kilo,quantity_units=nothing,clip=true,atol=nothing,kwargs...)
    requests=_plot_requests(first(results),_plot_ydata(selection,ydata,()))
    normalized=Grammar.observation_requests(first(results),requests;complete_pairs=true)
    observed=observables(results,requests;complete_pairs=true,frequency_unit=freq_unit,
        length_unit,quantity_units,clip,atol)
    observed_reference=reference===nothing || reference isa Grammar.ObservedResult ? reference :
        Grammar.ObservedResult(reference,requests;complete_pairs=true,frequency_unit=freq_unit,
            length_unit,quantity_units,clip,atol)
    return plot(observed;ydata=normalized.displayed,reference=observed_reference,problem,formulations,kwargs...)
end

function plot(results::LineCableModels.AbstractUncertaintyResult,selection=nothing;
        ydata=nothing,point=nothing,reference=nothing,freq_unit=:base,length_unit=:kilo,
        quantity_units=nothing,clip=true,atol=nothing,kwargs...)
    requests=_plot_requests(results,_plot_ydata(selection,ydata,()))
    observed=observables(results,requests;complete_pairs=true,frequency_unit=freq_unit,
        length_unit,quantity_units,clip,atol)
    retained=point===nothing ? observed : observed[point isa Integer ? [point] : point]
    observed_reference=reference===nothing || reference isa Grammar.ObservedResult ? reference :
        Grammar.ObservedResult(reference,1,requests;complete_pairs=true,frequency_unit=freq_unit,
            length_unit,quantity_units,clip,atol)
    return plot(retained;ydata=requests,reference=observed_reference,kwargs...)
end

"""Render a report's retained candidates and separate observed reference."""
function plot(artifact::ReportArtifact,selection=nothing;ydata=nothing,title_prefix=nothing,kwargs...)
    point=artifact.observed isa Grammar.ObservedResult ? artifact.observed : first(artifact.observed)
    prefix=title_prefix===nothing ? get(get(point.timings,:context,(;)),:id,nothing) : title_prefix
    return plot(artifact.observed,selection;ydata,reference=artifact.reference,title_prefix=prefix,kwargs...)
end
Makie.plot(artifact::ReportArtifact,args...;kwargs...) = plot(artifact,args...;kwargs...)
Makie.plot(results::LineCableModels.AbstractResultSpace,args...;kwargs...) = plot(results,args...;kwargs...)
