module PlotBuilder

using Base: @kwdef

import ..UnitHandler: Units, QuantityTag, get_label, get_symbol, display_unit, scale_factor

export PlotAxis, PlotRenderer,
	plot_kind, enable_logscale, dispatch_on,
	axis_quantity, axis_unit, axis_label,
	default_figsize,
	parse_kwargs, resolve_input, build_payloads, build_renderer

# Submodule `BackendHandler`
include("backendhandler/BackendHandler.jl")
using .BackendHandler

# Submodule `PlotUIComponents`
include("plotuicomponents/PlotUIComponents.jl")
using .PlotUIComponents


include("types.jl")
include("traits.jl")
include("plotaxis.jl")
include("parse.jl")
include("plotseries.jl")
include("plotarea.jl")
include("plotfigure.jl")


include("plotspecs.jl")


"""
	build_renderer(::Type{S}, obj; kwargs...) where {S<:AbstractPlotSpec}

High-level API: from domain object + keyword arguments to a PlotRenderer.

Checks that the object type is compatible with `dispatch_on(S)` and then
runs:

	parse_kwargs(S, obj; kwargs...)  → raw
	resolve_input(S, raw)            → nt
	build_figures(S, nt)             → figs

`build_figures` returns a vector of PlotFigure values; `build_renderer`
wraps them into a PlotRenderer that the UI layer will later assemble into
actual windows/layouts.
"""
function build_renderer(::Type{S}, obj; kwargs...) where {S <: AbstractPlotSpec}
	Tdispatch = dispatch_on(S)
	obj isa Tdispatch ||
		Base.error("Spec $(S) cannot dispatch on $(typeof(obj)); expected $(Tdispatch)")

	raw  = parse_kwargs(S, obj; kwargs...)
	norm = resolve_input(S, raw)
	figs = build_figures(S, norm)  # ::Vector{PlotFigure}

	return PlotRenderer(S, figs)
end


# """
# Final construction step: from normalized NamedTuple to concrete payloads.

# Return a `Vector{<:NamedTuple}` where each NamedTuple is the full Makie
# payload for a *single* plot primitive, with keys as close as possible
# to Makie’s own vocabulary, e.g.:

# 	(
# 		xdata   = ::AbstractVector,
# 		ydata   = ::AbstractVector,
# 		zdata   = ::Union{AbstractArray,Nothing},
# 		xaxis   = ::PlotAxis,
# 		yaxis   = ::PlotAxis,
# 		zaxis   = ::Union{PlotAxis,Nothing},
# 		title   = ::String,
# 		legend  = ::Vector{String},
# 		kwargs  = ::NamedTuple,   # color, linestyle, markersize, ...
# 	)

# Specs may override this when they need multiple primitives, grids, etc.
# """
# # --------------------------------------------------------------------------
# # Generic payload builder
# # --------------------------------------------------------------------------

# function build_payloads(::Type{S}, nt::NamedTuple) where {S <: AbstractPlotSpec}
# 	dims    = geom_axes(S)
# 	dims    = dims isa Tuple ? dims : (dims,)
# 	backend = nt.backend

# 	axes  = build_axes(S, nt)
# 	xaxis = axes.xaxis
# 	yaxis = axes.yaxis
# 	zaxis = axes.zaxis

# 	xdata = nothing
# 	ydata = nothing
# 	zdata = nothing

# 	:x in dims && (xdata = axis_data(S, :x, nt, xaxis))
# 	:y in dims && (ydata = axis_data(S, :y, nt, yaxis))
# 	:z in dims && (zdata = axis_data(S, :z, nt, zaxis))

# 	title  = default_title(S, nt)
# 	legend = legend_labels(S, nt)

# 	payload = (
# 		xdata  = xdata,
# 		ydata  = ydata,
# 		zdata  = zdata,
# 		xaxis  = xaxis,
# 		yaxis  = yaxis,
# 		zaxis  = zaxis,
# 		title  = title,
# 		legend = legend,
# 		kwargs = backend,
# 	)

# 	return [payload]
# end

# """
# High-level API: from domain object + kwargs to a vector of PlotRenderer.

# Checks that the object type is compatible with `dispatch_on(S)` and then
# runs:

# 	parse_kwargs → resolve_input → build_payloads

# `build_payloads` returns vector of payload NamedTuples; `build_renderer` wraps them into
# `PlotRenderer` values that the UI layer will later assemble into actual
# windows/layouts.
# """
# function build_renderer(::Type{S}, obj; kwargs...) where {S <: AbstractPlotSpec}
# 	Tdispatch = dispatch_on(S)
# 	obj isa Tdispatch ||
# 		Base.error("Spec $(S) cannot dispatch on $(typeof(obj)); expected $(Tdispatch)")

# 	raw      = parse_kwargs(S, obj; kwargs...)
# 	norm     = resolve_input(S, raw)
# 	payloads = build_payloads(S, norm)  # ::Vector{<:NamedTuple}

# 	renderers = Vector{PlotRenderer}(undef, length(payloads))
# 	@inbounds for (k, p) in pairs(payloads)
# 		renderers[k] = PlotRenderer(S, p)
# 	end
# 	return renderers
# end

end # module PlotBuilder
