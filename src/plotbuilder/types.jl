abstract type AbstractPlotSpec end

"""
Axis is the fully decided axis descriptor used by plot areas (panels).

- `dim`      : axis dimension (e.g. :x, :y, :z)
- `quantity` : semantic quantity tag (from UnitHandler)
- `units`    : unit system for this axis
- `label`    : full label text, including unit symbol
- `scale`    : :linear or :log10
"""
struct Axis
	dim::Symbol
	quantity::QuantityTag
	units::Units
	label::String
	scale::Symbol
end

# --------------------------------------------------------------------------
# Payload hierarchy: series → area → figure → renderer
# --------------------------------------------------------------------------
"""
	Dataseries

Single plot primitive (one Makie call).

Fields:
- `kind`   : plotting primitive kind (:line, :scatter, :hist, :heatmap, ...)
- `xdata`  : x values \\[dimensionless or scaled\\]
- `ydata`  : y values \\[dimensionless or scaled\\]
- `zdata`  : z values if applicable, otherwise `nothing`
- `label`  : legend entry for this series, or `nothing` for no legend
"""
struct Dataseries
	kind  :: Symbol
	xdata :: Union{Nothing, AbstractVector{<:Number}}
	ydata :: Union{Nothing, AbstractArray{<:Number}}
	zdata :: Union{Nothing, AbstractArray{<:Number}}
	label :: Union{Nothing, String}
end

"""
	Panel

One plot panel / axis system.

All Dataseries inside a Panel share the same x/y(/z) axes. The `key`
field encodes the semantic facet this area represents (e.g. \\(i,j\\),
quantity, frequency).

Fields:
- `xaxis`  : x Axis or `nothing` if unused
- `yaxis`  : y Axis or `nothing` if unused
- `zaxis`  : z Axis or `nothing` if unused
- `title`  : panel title
- `series` : vector of Dataseries
- `key`    : NamedTuple identifying the facet, or empty `NamedTuple` if none
"""
struct Panel
	xaxis  :: Union{Nothing, Axis}
	yaxis  :: Union{Nothing, Axis}
	zaxis  :: Union{Nothing, Axis}
	title  :: String
	series :: Vector{Dataseries}
	key    :: NamedTuple
end

"""
	Figure

One logical figure / window.

Fields:
- `title`      : optional figure-level title (may be empty)
- `size`       : (width, height) in pixels
- `layout`     : layout spec (:windows, :grid, :tabbed, ...)
- `panels` : vector of Panel values contained in this figure
- `kwargs`     : figure-level backend options (e.g. figsize)
"""
struct Figure
	title  :: String
	size   :: Tuple{Int, Int}
	layout :: Symbol
	panels :: Vector{Panel}
	kwargs :: NamedTuple
end

"""
	PlotRenderer{S}

Final product of the grammar pipeline for spec type `S`.

Carries only:
- the spec type (the grammar definition),
- a vector of Figure payloads.

The rendering backend (Makie) should only see PlotRenderer values and must
never touch domain objects or grammar logic.
"""
struct PlotRenderer{S <: AbstractPlotSpec}
	spec_type :: Type{S}
	figures   :: Vector{Figure}
end
