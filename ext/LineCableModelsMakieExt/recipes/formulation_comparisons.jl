import LineCableModels.ReportBuilder: ReportArtifact

"""Construct observations, then overlay selected original study points."""
function plot(results::LineCableModels.ParametricResult,selection=nothing;
        problem=nothing,formulations=nothing,reference=nothing,kwargs...)
    acquisition,presentation=_plot_observation_options(kwargs)
    requests=Commons.observation_selection(first(results),selection)
    normalized=Commons.observation_requests(first(results),requests;complete_pairs=true)
    observed=observables(results,requests;complete_pairs=true,acquisition...)
    validate(reference,plot)
    if reference!==nothing && !(reference isa Commons.ObservedResult)
        reference=Commons.ObservedResult(reference,requests;complete_pairs=true,acquisition...)
    end
    return plot(observed,normalized.displayed;reference,problem,formulations,presentation...)
end

function plot(results::LineCableModels.AbstractUncertaintyResult,selection=nothing;
        point=nothing,reference=nothing,kwargs...)
    acquisition,presentation=_plot_observation_options(kwargs)
    requests=Commons.observation_selection(results,selection)
    observed=observables(results,requests;complete_pairs=true,acquisition...)
    # A number selects one point, as a collection of one.
    retained=point===nothing ? observed : observed[[index for index in point]]
    validate(reference,plot)
    if reference!==nothing && !(reference isa Commons.ObservedResult)
        reference=Commons.ObservedResult(reference,requests;complete_pairs=true,acquisition...)
    end
    return plot(retained,requests;reference,presentation...)
end

"""Render a report's retained observations and separate observed reference."""
function plot(artifact::ReportArtifact,selection=nothing;title_prefix=nothing,reference=artifact.reference,kwargs...)
    point=artifact.observed isa Commons.ObservedResult ? artifact.observed : first(artifact.observed)
    prefix=title_prefix===nothing ? get(get(point.timings,:context,(;)),:id,nothing) : title_prefix
    return plot(artifact.observed,selection;reference,title_prefix=prefix,kwargs...)
end
Makie.plot(artifact::ReportArtifact,args...;kwargs...) = plot(artifact,args...;kwargs...)
Makie.plot(results::LineCableModels.AbstractResultSpace,args...;kwargs...) = plot(results,args...;kwargs...)
