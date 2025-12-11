module PlotBuilder

using Base: @kwdef

import ..UnitHandler: Units, QuantityTag, get_label, get_symbol, display_unit, scale_factor

export build_renderer, PlotRenderer

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

end # module PlotBuilder
