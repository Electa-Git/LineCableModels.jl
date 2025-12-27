module PlotBuilder

using Base: @kwdef

import ..UnitHandler: Units, QuantityTag, get_label, get_symbol, display_unit, scale_factor
import ..Commons: PhaseDomain, ModalDomain, domain

export make_render, RenderSpec

# Submodule `BackendHandler`
include("backendhandler/BackendHandler.jl")
using .BackendHandler

# Submodule `PlotUIComponents` - WILL BE DEPRECATED SOON
include("plotuicomponents/PlotUIComponents.jl")
using .PlotUIComponents

include("types.jl")
include("traits.jl")
include("axisspec.jl")
include("parse.jl")
include("seriesspec.jl")
include("viewspec.jl")
include("pagespec.jl")

# Submodule `UIComponents`
include("uicomponents/UIComponents.jl")

include("plotspecs.jl")



"""
	make_render(::Type{S}, obj; kwargs...) where {S<:AbstractPlotSpec}

High-level API: from domain object + keyword arguments to a RenderSpec.

Checks that the object type is compatible with `dispatch_on(S)` and then
runs:

	parse_kwargs(S, obj; kwargs...)  → raw
	resolve_input(S, raw)            → nt
	make_pages(S, nt)             → figs

`make_pages` returns a vector of PageSpec values; `make_render`
wraps them into a RenderSpec that the UI layer will later assemble into
actual windows/layouts.
"""

function make_render(::Type{S}, obj; kwargs...) where {S <: AbstractPlotSpec}
	Tdispatch = dispatch_on(S)
	obj isa Tdispatch ||
		Base.error("Spec $(S) cannot dispatch on $(typeof(obj)); expected $(Tdispatch)")

	raw  = parse_kwargs(S, obj; kwargs...)
	norm = resolve_input(S, raw)
	pags = make_pages(S, norm)  # ::Vector{PageSpec}

	return RenderSpec(S, pags)
end

end # module PlotBuilder
