"""
    LineCableModelsMakieExt

Build Makie figures and controls from LineCableModels plot recipes.
"""
module LineCableModelsMakieExt

using LineCableModels
using Makie

const PlotBuilder = LineCableModels.PlotBuilder

function current_backend_symbol()
    backend = Makie.current_backend()
    backend isa Module || return :none
    name = nameof(backend)
    name === :CairoMakie && return :cairo
    name === :GLMakie && return :gl
    name === :WGLMakie && return :wgl
    return :unknown
end

renderfig(figure) = display(figure)

include("UIComponents.jl")
using .UIComponents

import LineCableModels.Engine: plot
import LineCableModels.DataModel: preview, show_material_scale

function _scale_symbol(value)
    value isa Symbol && return value
    value === Makie.identity && return :linear
    value === Makie.log10 && return :log10
    throw(ArgumentError("axis scale must be :linear, :log10, Makie.identity, or Makie.log10"))
end

function plot(
        object::LineCableModels.SeriesImpedance,
        frequencies,
        quantities::Tuple = ();
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        xscale = :linear,
        yscale = :linear,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotDefinition,
        object;
        frequencies,
        quantities,
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        kwargs...
    )
    return UIComponents.build(render_spec; backend, display = display_plot, controls)
end

function Makie.plot(
        object::LineCableModels.SeriesImpedance,
        frequencies,
        quantities::Tuple = ();
        kwargs...
)
    plot(object, frequencies, quantities; kwargs...)
end

function plot(
        object::LineCableModels.ShuntAdmittance,
        frequencies,
        quantities::Tuple = ();
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        xscale = :linear,
        yscale = :linear,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotDefinition,
        object;
        frequencies,
        quantities,
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        kwargs...
    )
    return UIComponents.build(render_spec; backend, display = display_plot, controls)
end

function Makie.plot(
        object::LineCableModels.ShuntAdmittance,
        frequencies,
        quantities::Tuple = ();
        kwargs...
)
    plot(object, frequencies, quantities; kwargs...)
end

function plot(
        parameters::LineCableModels.LineParameters,
        quantities::Tuple = ();
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        xscale = :linear,
        yscale = :linear,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.Engine.LineParameterPlotDefinition,
        parameters;
        quantities,
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        kwargs...
    )
    return UIComponents.build(render_spec; backend, display = display_plot, controls)
end

function Makie.plot(
        parameters::LineCableModels.LineParameters,
        quantities::Tuple = ();
        kwargs...
)
    plot(parameters, quantities; kwargs...)
end

function plot(
        first::LineCableModels.LineParameters,
        second::LineCableModels.LineParameters,
        rest::LineCableModels.LineParameters...;
        legend,
        quantities::Tuple = (),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        xscale = :linear,
        yscale = :linear,
        kwargs...
)
    parameters = (first, second, rest...)
    render_spec = PlotBuilder.make_render(
        LineCableModels.Engine.LineParametersBenchmarkPlotDefinition,
        parameters;
        legend,
        quantities,
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        kwargs...
    )
    return UIComponents.build(render_spec; backend, display = display_plot, controls)
end

function Makie.plot(
        first::LineCableModels.LineParameters,
        second::LineCableModels.LineParameters,
        rest::LineCableModels.LineParameters...;
        kwargs...
)
    return plot(first, second, rest...; kwargs...)
end

function plot(
        result::LineCableModels.MonteCarloResult,
        selector::Function = LineCableModels.R;
        ijk = nothing,
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.UQ.MCDistributionPlotDefinition,
        result;
        selector,
        ijk,
        kwargs...
    )
    return only(UIComponents.build(render_spec; backend, display = display_plot, controls))
end

function Makie.plot(
        result::LineCableModels.MonteCarloResult,
        selector::Function = LineCableModels.R;
        kwargs...
)
    return plot(result, selector; kwargs...)
end

function preview(
        design::LineCableModels.DataModel.CableDesign;
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.DataModel.CablePreviewPlotDefinition,
        design;
        kwargs...
    )
    return only(UIComponents.build(render_spec; backend, display = display_plot, controls))
end

function preview(
        system::LineCableModels.DataModel.LineCableSystem;
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.DataModel.SystemPreviewPlotDefinition,
        system;
        kwargs...
    )
    return only(UIComponents.build(render_spec; backend, display = display_plot, controls))
end

function show_material_scale(
        ; backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    render_spec = PlotBuilder.make_render(
        LineCableModels.DataModel.MaterialScalePlotDefinition,
        nothing;
        kwargs...
    )
    return only(UIComponents.build(render_spec; backend, display = display_plot, controls))
end

end # module LineCableModelsMakieExt
