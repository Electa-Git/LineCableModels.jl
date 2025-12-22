abstract type AbstractPlotSpec end

"""
PBAxis is the fully decided axis descriptor used by plot areas (views).

- `dim`      : axis dimension (e.g. :x, :y, :z)
- `quantity` : semantic quantity tag (from UnitHandler)
- `units`    : unit system for this axis
- `label`    : full label text, including unit symbol
- `scale`    : :linear or :log10
"""
struct PBAxis
	dim::Symbol
	quantity::QuantityTag
	units::Units
	label::String
	scale::Symbol
end

# --------------------------------------------------------------------------
# Payload hierarchy: series → view → figure → renderer
# --------------------------------------------------------------------------
"""
	PBSeries

Single plot primitive (one Makie call).

Fields:
- `kind`   : plotting primitive kind (:line, :scatter, :hist, :heatmap, ...)
- `xdata`  : x values \\[dimensionless or scaled\\]
- `ydata`  : y values \\[dimensionless or scaled\\]
- `zdata`  : z values if applicable, otherwise `nothing`
- `label`  : legend entry for this series, or `nothing` for no legend
"""
struct PBSeries
	kind  :: Symbol
	xdata :: Union{Nothing, AbstractVector{<:Number}}
	ydata :: Union{Nothing, AbstractArray{<:Number}}
	zdata :: Union{Nothing, AbstractArray{<:Number}}
	label :: Union{Nothing, String}
end

"""
	PBView

One plot view / axis system.

All PBSeries inside a PBView share the same x/y(/z) axes. The `key`
field encodes the semantic facet this area represents (e.g. \\(i,j\\),
quantity, frequency).

Fields:
- `xaxis`  : x PBAxis or `nothing` if unused
- `yaxis`  : y PBAxis or `nothing` if unused
- `zaxis`  : z PBAxis or `nothing` if unused
- `title`  : view title
- `series` : vector of PBSeries
- `key`    : NamedTuple identifying the facet, or empty `NamedTuple` if none
"""
struct PBView
	xaxis  :: Union{Nothing, PBAxis}
	yaxis  :: Union{Nothing, PBAxis}
	zaxis  :: Union{Nothing, PBAxis}
	title  :: String
	series :: Vector{PBSeries}
	key    :: NamedTuple
end

"""
	PBFigure

One logical figure / window.

Fields:
- `title`      : optional figure-level title (may be empty)
- `size`       : (width, height) in pixels
- `layout`     : layout spec (:windows, :grid, :tabbed, ...)
- `views`      : vector of PBView values contained in this figure
- `kwargs`     : figure-level backend options (e.g. figsize)
"""
struct PBFigure
	title  :: String
	size   :: Tuple{Int, Int}
	layout :: Symbol
	views  :: Vector{PBView}
	kwargs :: NamedTuple
end

"""
	PBRenderer{S}

Final product of the grammar pipeline for spec type `S`.

Carries only:
- the spec type (the grammar definition),
- a vector of PBFigure payloads.

The rendering backend (Makie) should only see PBRenderer values and must
never touch domain objects or grammar logic.
"""
struct PBRenderer{S <: AbstractPlotSpec}
	spec    :: Type{S}
	figures :: Vector{PBFigure}
end
