import LineCableModels.ReportBuilder: ReportArtifact

"""Construct observations, then overlay selected original study points."""
function plot(results::LineCableModels.ParametricResult,selection=nothing;
        ydata=nothing,problem=nothing,formulations=nothing,reference=nothing,kwargs...)
    acquisition,presentation=_plot_observation_options(kwargs)
    requests=Grammar.observation_selection(first(results),_plot_ydata(selection,ydata,()))
    normalized=Grammar.observation_requests(first(results),requests;complete_pairs=true)
    observed=observables(results,requests;complete_pairs=true,acquisition...)
    if reference!==nothing && !(reference isa Grammar.ObservedResult)
        reference isa _PrimarySource || throw(ArgumentError("construct an atomic ObservedResult for a reference collection"))
        reference=Grammar.ObservedResult(reference,requests;complete_pairs=true,acquisition...)
    end
    return plot(observed;ydata=normalized.displayed,reference,problem,formulations,presentation...)
end

function plot(results::LineCableModels.AbstractUncertaintyResult,selection=nothing;
        ydata=nothing,point=nothing,reference=nothing,kwargs...)
    acquisition,presentation=_plot_observation_options(kwargs)
    requests=Grammar.observation_selection(results,_plot_ydata(selection,ydata,()))
    observed=observables(results,requests;complete_pairs=true,acquisition...)
    retained=point===nothing ? observed : observed[point isa Integer ? [point] : point]
    if reference!==nothing && !(reference isa Grammar.ObservedResult)
        reference isa _PrimarySource || throw(ArgumentError("construct an atomic ObservedResult for a reference collection"))
        reference=Grammar.ObservedResult(reference,requests;complete_pairs=true,acquisition...)
    end
    return plot(retained;ydata=requests,reference,presentation...)
end

"""Render a report's retained candidates and separate observed reference."""
function plot(artifact::ReportArtifact,selection=nothing;ydata=nothing,title_prefix=nothing,reference=artifact.reference,kwargs...)
    point=artifact.observed isa Grammar.ObservedResult ? artifact.observed : first(artifact.observed)
    prefix=title_prefix===nothing ? get(get(point.timings,:context,(;)),:id,nothing) : title_prefix
    return plot(artifact.observed,selection;ydata,reference,title_prefix=prefix,kwargs...)
end
Makie.plot(artifact::ReportArtifact,args...;kwargs...) = plot(artifact,args...;kwargs...)
Makie.plot(results::LineCableModels.AbstractResultSpace,args...;kwargs...) = plot(results,args...;kwargs...)
