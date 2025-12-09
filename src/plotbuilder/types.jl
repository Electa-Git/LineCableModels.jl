abstract type AbstractPlotSpec end

"""
PlotAxis is the fully decided axis descriptor used by plot areas.

- `dim`      : axis dimension (e.g. :x, :y, :z)
- `quantity` : semantic quantity tag (from UnitHandler)
- `units`    : unit system for this axis
- `label`    : full label text, including unit symbol
- `scale`    : :linear or :log10
"""
struct PlotAxis
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
	PlotSeries

Single plot primitive (one Makie call).

Fields:
- `kind`   : plotting primitive kind (:line, :scatter, :hist, :heatmap, ...)
- `xdata`  : x values \\[dimensionless or scaled\\]
- `ydata`  : y values \\[dimensionless or scaled\\]
- `zdata`  : z values if applicable, otherwise `nothing`
- `label`  : legend entry for this series, or `nothing` for no legend
- `kwargs` : backend kwargs, passed as-is to the rendering backend
"""
struct PlotSeries
	kind   :: Symbol
	xdata  :: Union{Nothing, AbstractVector{<:Number}}
	ydata  :: Union{Nothing, AbstractVector{<:Number}}
	zdata  :: Union{Nothing, AbstractVector{<:Number}}
	label  :: Union{Nothing, String}
	kwargs :: NamedTuple
end

"""
	PlotArea

One plot panel / axis system.

All PlotSeries inside a PlotArea share the same x/y(/z) axes. The `key`
field encodes the semantic facet this area represents (e.g. \\(i,j\\),
quantity, frequency).

Fields:
- `xaxis`  : x PlotAxis or `nothing` if unused
- `yaxis`  : y PlotAxis or `nothing` if unused
- `zaxis`  : z PlotAxis or `nothing` if unused
- `title`  : panel title
- `series` : vector of PlotSeries
- `key`    : NamedTuple identifying the facet, or empty `NamedTuple` if none
"""
struct PlotArea
	xaxis  :: Union{Nothing, PlotAxis}
	yaxis  :: Union{Nothing, PlotAxis}
	zaxis  :: Union{Nothing, PlotAxis}
	title  :: String
	series :: Vector{PlotSeries}
	key    :: NamedTuple
end

"""
	PlotFigure

One logical figure / window.

Fields:
- `title`      : optional figure-level title (may be empty)
- `layout`     : layout hint (:single, :grid, :free, ...)
- `plot_areas` : vector of PlotArea values contained in this figure
- `kwargs`     : figure-level backend options (e.g. figsize)
"""
struct PlotFigure
	title      :: String
	layout     :: Symbol
	plot_areas :: Vector{PlotArea}
	kwargs     :: NamedTuple
end

"""
	PlotRenderer{S}

Final product of the grammar pipeline for spec type `S`.

Carries only:
- the spec type (the grammar definition),
- a vector of PlotFigure payloads.

The rendering backend (Makie) should only see PlotRenderer values and must
never touch domain objects or grammar logic.
"""
struct PlotRenderer{S <: AbstractPlotSpec}
	spec_type :: Type{S}
	figures   :: Vector{PlotFigure}
end
