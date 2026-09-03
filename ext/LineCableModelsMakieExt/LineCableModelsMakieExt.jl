"""
    LineCableModelsMakieExt

Add compact high-level LineCableModels plotting methods to native Makie.
"""
module LineCableModelsMakieExt

using LineCableModels
using Makie
using LinearAlgebra: diag
using Printf: @sprintf
using Colors: HSL, HSV, RGB, RGBA, alpha, blue, green, red
using Dates
using Statistics: mean

const Units = LineCableModels.Units
import LineCableModels.Grammar:
                                observation_indices, observation_request, request_identity,
                                request_indices,
                                request_quantity, unit_targets, validate_observables

function current_backend_symbol()
    backend = Makie.current_backend()
    backend isa Module || return :none
    name = nameof(backend)
    name === :CairoMakie && return :cairo
    name === :GLMakie && return :gl
    name === :WGLMakie && return :wgl
    return :unknown
end

include("recipes/line_data.jl")
include("recipes/comparison_data.jl")
include("material_colors.jl")
include("recipes/preview_types.jl")
include("recipes/preview_data.jl")
include("shell.jl")
include("recipes/line_facets.jl")
include("recipes/preview_render.jl")
include("recipes/publication_render.jl")
include("montecarlo.jl")
include("native_export.jl")

import LineCableModels.PlotBuilder: plot, preview, show_material_scale

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

function _is_line_index_selector(value)
    value isa Integer || value isa AbstractRange ||
        value isa AbstractVector{<:Integer} ||
        value isa Colon
end

function _is_line_observation_request(source, value)
    value isa Tuple || return false
    resolved = try
        LineCableModels.Grammar.observation_request(source, value)
    catch error
        error isa ArgumentError || rethrow()
        return false
    end
    expected = resolved.identity isa Tuple && last(resolved.identity) === diag ? 2 : 3
    return length(resolved.indices) == expected &&
           all(_is_line_index_selector, resolved.indices)
end

function _expand_line_observation_request(source, request)
    resolved = LineCableModels.Grammar.observation_request(source, request)
    identity = resolved.identity
    indices = resolved.indices
    if identity === LineCableModels.Z
        return ((LineCableModels.R, indices...), (LineCableModels.X, indices...))
    elseif identity === LineCableModels.Y
        return ((LineCableModels.G, indices...), (LineCableModels.B, indices...))
    elseif identity == (LineCableModels.Z, diag)
        return (
            (LineCableModels.R, diag, indices...),
            (LineCableModels.X, diag, indices...)
        )
    elseif identity == (LineCableModels.Y, diag)
        return (
            (LineCableModels.G, diag, indices...),
            (LineCableModels.B, diag, indices...)
        )
    end
    return (request,)
end

_full_line_request(selector::Function) = (selector, Colon(), Colon(), Colon())
function _full_line_request(selector::Function, transform::Function)
    (selector, transform, Colon(), Colon(), Colon())
end
_full_diagonal_request(selector::Function) = (selector, diag, Colon(), Colon())

function _domain_line_request(source, selector::Function)
    LineCableModels.domain(source) === LineCableModels.ModalDomain &&
        return _full_diagonal_request(selector)
    return _full_line_request(selector)
end

function _line_selector_requests(source::LineCableModels.LineParameters, selector::Function)
    selector === LineCableModels.Z &&
        return (_domain_line_request(source, LineCableModels.R),
            _domain_line_request(source, LineCableModels.X))
    selector === LineCableModels.Y &&
        return (_domain_line_request(source, LineCableModels.G),
            _domain_line_request(source, LineCableModels.B))
    selector === real && return (_domain_line_request(source, LineCableModels.R),
        _domain_line_request(source, LineCableModels.G))
    selector === imag && return (_domain_line_request(source, LineCableModels.X),
        _domain_line_request(source, LineCableModels.B))
    selector === abs && return (_full_line_request(LineCableModels.Z, abs),
        _full_line_request(LineCableModels.Y, abs))
    selector === angle && return (_full_line_request(LineCableModels.Z, angle),
        _full_line_request(LineCableModels.Y, angle))
    return (_domain_line_request(source, selector),)
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

function _default_line_selection(::LineCableModels.LineParameters)
    (
        LineCableModels.R, LineCableModels.X, LineCableModels.G, LineCableModels.B)
end
function _default_line_selection(::LineCableModels.SeriesImpedance)
    (LineCableModels.R, LineCableModels.X)
end
function _default_line_selection(::LineCableModels.ShuntAdmittance)
    (LineCableModels.G, LineCableModels.B)
end

function _line_plot_requests(source::_LineSource, selection)
    selected = selection === nothing || selection == () ?
               _default_line_selection(source) : selection
    _is_line_observation_request(source, selected) &&
        return _expand_line_observation_request(source, selected)
    selected isa Function && return _line_selector_requests(source, selected)
    selected isa Tuple || throw(ArgumentError(
        "line selection must be a selector, observable request, or tuple of them",
    ))
    return Tuple(request
    for item in selected
    for request in (_is_line_observation_request(source, item) ?
         _expand_line_observation_request(source, item) :
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
    return _addon_semantic_line_plots(
        object;
        frequencies,
        requests = _line_plot_requests(object, selection),
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        backend,
        display_plot,
        controls,
        kwargs...
    )
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
    return _addon_semantic_line_plots(
        object;
        frequencies,
        requests = _line_plot_requests(object, selection),
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        backend,
        display_plot,
        controls,
        kwargs...
    )
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
    return _addon_semantic_line_plots(
        parameters;
        requests = _line_plot_requests(parameters, selection),
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        backend,
        display_plot,
        controls,
        kwargs...
    )
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
        rest...;
        series_labels = nothing,
        requests = (),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        xscale = :linear,
        yscale = :linear,
        kwargs...
)
    sources = LineCableModels.LineParameters[first, second]
    trailing_selection = nothing
    for item in rest
        if item isa LineCableModels.LineParameters && trailing_selection === nothing
            push!(sources, item)
        elseif trailing_selection === nothing
            trailing_selection = item
        else
            throw(ArgumentError(
                "a line comparison accepts sources followed by at most one observation selection",
            ))
        end
    end
    requests == () || trailing_selection === nothing ||
        throw(ArgumentError(
            "use either a trailing observation selection or the requests keyword, not both",
        ))
    labels = series_labels
    labels === nothing && (labels = Tuple("Result $index" for index in eachindex(sources)))
    selection = trailing_selection === nothing ? requests : trailing_selection
    normalized = _line_plot_requests(first, selection)
    return _addon_line_pages(
        Tuple(sources);
        series_labels = labels,
        requests = normalized,
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

function plot(
        sources::NamedTuple,
        selection = ();
        series_labels = nothing,
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
    labels = if series_labels !== nothing
        series_labels
    else
        Tuple(String(key) for key in keys(sources))
    end
    return _addon_line_pages(
        parameters;
        series_labels = labels,
        requests = _line_plot_requests(first(parameters), selection),
        xscale = _scale_symbol(xscale),
        yscale = _scale_symbol(yscale),
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

function Makie.plot(sources::NamedTuple, selection = (); kwargs...)
    plot(sources, selection; kwargs...)
end

function Makie.plot(
        first::LineCableModels.LineParameters,
        second::LineCableModels.LineParameters,
        rest...;
        kwargs...
)
    return plot(first, second, rest...; kwargs...)
end

function preview(
        design::LineCableModels.DataModel.CableDesign;
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    return _addon_preview(
        design;
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

# The public extension method consumes backend/display choices and forwards the
# remaining preview options unchanged. DataModel retains only detached geometry
# and material attributes; Makie objects and backend state stay here.
function preview(
        designs::AbstractVector{<:LineCableModels.DataModel.CableDesign};
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    return _addon_preview(
        designs;
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

function preview(
        system::LineCableModels.DataModel.LineCableSystem;
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    return _addon_preview(
        system;
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

function show_material_scale(
        ; backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    return _addon_material_scale(;
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

function plot(
        publication::LineCableModels.Grammar.ObservationPublication;
        kwargs...
)
    return _addon_publication_plot(publication; kwargs...)
end

function Makie.plot(publication::LineCableModels.Grammar.ObservationPublication; kwargs...)
    plot(publication; kwargs...)
end

end # module LineCableModelsMakieExt
