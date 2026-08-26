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

const _LineSource = Union{
    LineCableModels.LineParameters,
    LineCableModels.SeriesImpedance,
    LineCableModels.ShuntAdmittance
}

_is_line_index_selector(value) = value isa Integer || value isa AbstractRange ||
                                 value isa AbstractVector{<:Integer} ||
                                 value isa Colon

function _is_line_observation_request(value)
    value isa Tuple || return false
    length(value) in (4, 5) || return false
    value[1] isa Function || return false
    offset = length(value) == 4 ? 1 : 2
    length(value) == 5 && !(value[2] isa Function) && return false
    return all(_is_line_index_selector, value[(offset + 1):end])
end

_full_line_request(selector::Function) = (selector, Colon(), Colon(), Colon())
_full_line_request(selector::Function, transform::Function) =
    (selector, transform, Colon(), Colon(), Colon())

function _line_selector_requests(::LineCableModels.LineParameters, selector::Function)
    selector === LineCableModels.Z && return (_full_line_request(LineCableModels.R),
        _full_line_request(LineCableModels.X))
    selector === LineCableModels.Y && return (_full_line_request(LineCableModels.G),
        _full_line_request(LineCableModels.B))
    selector === real && return (_full_line_request(LineCableModels.R),
        _full_line_request(LineCableModels.G))
    selector === imag && return (_full_line_request(LineCableModels.X),
        _full_line_request(LineCableModels.B))
    selector === abs && return (_full_line_request(LineCableModels.Z, abs),
        _full_line_request(LineCableModels.Y, abs))
    selector === angle && return (_full_line_request(LineCableModels.Z, angle),
        _full_line_request(LineCableModels.Y, angle))
    return (_full_line_request(selector),)
end

function _line_selector_requests(::LineCableModels.SeriesImpedance, selector::Function)
    selector === LineCableModels.Z && return (_full_line_request(LineCableModels.R),
        _full_line_request(LineCableModels.X))
    selector === real && return (_full_line_request(LineCableModels.R),)
    selector === imag && return (_full_line_request(LineCableModels.X),)
    selector === abs && return (_full_line_request(LineCableModels.Z, abs),)
    selector === angle && return (_full_line_request(LineCableModels.Z, angle),)
    return (_full_line_request(selector),)
end

function _line_selector_requests(::LineCableModels.ShuntAdmittance, selector::Function)
    selector === LineCableModels.Y && return (_full_line_request(LineCableModels.G),
        _full_line_request(LineCableModels.B))
    selector === real && return (_full_line_request(LineCableModels.G),)
    selector === imag && return (_full_line_request(LineCableModels.B),)
    selector === abs && return (_full_line_request(LineCableModels.Y, abs),)
    selector === angle && return (_full_line_request(LineCableModels.Y, angle),)
    return (_full_line_request(selector),)
end

_default_line_selection(::LineCableModels.LineParameters) = (
    LineCableModels.R, LineCableModels.X, LineCableModels.G, LineCableModels.B)
_default_line_selection(::LineCableModels.SeriesImpedance) =
    (LineCableModels.R, LineCableModels.X)
_default_line_selection(::LineCableModels.ShuntAdmittance) =
    (LineCableModels.G, LineCableModels.B)

function _line_plot_requests(source::_LineSource, selection)
    selected = selection === nothing || selection == () ?
               _default_line_selection(source) : selection
    _is_line_observation_request(selected) && return (selected,)
    selected isa Function && return _line_selector_requests(source, selected)
    selected isa Tuple || throw(ArgumentError(
        "line selection must be a selector, observable request, or tuple of them",
    ))
    return Tuple(request
    for item in selected
    for request in (_is_line_observation_request(item) ? (item,) :
                    item isa Function ? _line_selector_requests(source, item) :
                    throw(ArgumentError("unsupported line selection $(repr(item))"))))
end

function plot(
        object::LineCableModels.SeriesImpedance,
        frequencies,
        selection = ();
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
        requests = _line_plot_requests(object, selection),
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        kwargs...
    )
    return UIComponents.build(render_spec; backend, display = display_plot, controls)
end

function Makie.plot(
        object::LineCableModels.SeriesImpedance,
        frequencies,
        selection = ();
        kwargs...
)
    plot(object, frequencies, selection; kwargs...)
end

function plot(
        object::LineCableModels.ShuntAdmittance,
        frequencies,
        selection = ();
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
        requests = _line_plot_requests(object, selection),
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        kwargs...
    )
    return UIComponents.build(render_spec; backend, display = display_plot, controls)
end

function Makie.plot(
        object::LineCableModels.ShuntAdmittance,
        frequencies,
        selection = ();
        kwargs...
)
    plot(object, frequencies, selection; kwargs...)
end

function plot(
        parameters::LineCableModels.LineParameters,
        selection = ();
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
        requests = _line_plot_requests(parameters, selection),
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        kwargs...
    )
    return UIComponents.build(render_spec; backend, display = display_plot, controls)
end

function Makie.plot(
        parameters::LineCableModels.LineParameters,
        selection = ();
        kwargs...
)
    plot(parameters, selection; kwargs...)
end

function plot(
        first::LineCableModels.LineParameters,
        second::LineCableModels.LineParameters,
        rest::LineCableModels.LineParameters...;
        legend,
        requests = (),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        xscale = :linear,
        yscale = :linear,
        kwargs...
)
    parameters = (first, second, rest...)
    normalized = _line_plot_requests(first, requests)
    render_spec = PlotBuilder.make_render(
        LineCableModels.Engine.LineParametersBenchmarkPlotDefinition,
        parameters;
        legend,
        requests = normalized,
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        kwargs...
    )
    return UIComponents.build(render_spec; backend, display = display_plot, controls)
end


function plot(
        sources::NamedTuple,
        selection = ();
        legend = nothing,
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        xscale = :linear,
        yscale = :linear,
        kwargs...
)
    parameters = Tuple(values(sources))
    length(parameters) >= 2 || throw(ArgumentError(
        "line comparison requires at least two named sources",
    ))
    all(parameter -> parameter isa LineCableModels.LineParameters, parameters) ||
        throw(ArgumentError("all comparison sources must be LineParameters"))
    labels = legend === nothing ? Tuple(String(key) for key in keys(sources)) : legend
    render_spec = PlotBuilder.make_render(
        LineCableModels.Engine.LineParametersBenchmarkPlotDefinition,
        parameters;
        legend = labels,
        requests = _line_plot_requests(first(parameters), selection),
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        kwargs...
    )
    return UIComponents.build(render_spec; backend, display = display_plot, controls)
end

Makie.plot(sources::NamedTuple, selection = (); kwargs...) =
    plot(sources, selection; kwargs...)

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
    # This extension method removes Makie-only execution keywords, then asks UQ
    # to compile the scientific request into one detached PlotBuilder page.
    # `selector` is positional for the familiar plotting call; `ijk` remains a
    # named coordinate because it applies only to matrix-valued line parameters.
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
