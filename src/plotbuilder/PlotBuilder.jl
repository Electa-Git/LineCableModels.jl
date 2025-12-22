module PlotBuilder

using Base: @kwdef

import ..UnitHandler: Units, QuantityTag, get_label, get_symbol, display_unit, scale_factor
import ..Commons: PhaseDomain, ModalDomain, domain

export build_renderer, PBRenderer

# Submodule `BackendHandler`
include("backendhandler/BackendHandler.jl")
using .BackendHandler

# Submodule `PlotUIComponents`
include("plotuicomponents/PlotUIComponents.jl")
using .PlotUIComponents


include("types.jl")
include("traits.jl")
include("axes.jl")
include("parse.jl")
include("dataseries.jl")
include("views.jl")
include("figures.jl")


include("plotspecs.jl")


"""
	build_renderer(::Type{S}, obj; kwargs...) where {S<:AbstractPlotSpec}

High-level API: from domain object + keyword arguments to a PBRenderer.

Checks that the object type is compatible with `dispatch_on(S)` and then
runs:

	parse_kwargs(S, obj; kwargs...)  → raw
	resolve_input(S, raw)            → nt
	build_figures(S, nt)             → figs

`build_figures` returns a vector of PBFigure values; `build_renderer`
wraps them into a PBRenderer that the UI layer will later assemble into
actual windows/layouts.
"""
function build_renderer(::Type{S}, obj; kwargs...) where {S <: AbstractPlotSpec}
	Tdispatch = dispatch_on(S)
	obj isa Tdispatch ||
		Base.error("Spec $(S) cannot dispatch on $(typeof(obj)); expected $(Tdispatch)")

	raw  = parse_kwargs(S, obj; kwargs...)
	norm = resolve_input(S, raw)
	figs = build_figures(S, norm)  # ::Vector{PBFigure}

	return PBRenderer(S, figs)
end

end # module PlotBuilder
