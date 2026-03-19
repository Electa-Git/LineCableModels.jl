abstract type AbstractPlotSpec end

"""
AxisSpec is the fully decided axis descriptor used by plot areas (views).

- `dim`      : axis dimension (e.g. :x, :y, :z)
- `quantity` : semantic quantity tag (from UnitHandler)
- `units`    : unit system for this axis
- `label`    : full label text, including unit symbol
- `scale`    : :linear or :log10
"""
struct AxisSpec
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
	SeriesSpec

Single plot primitive (one Makie call).

Fields:
- `kind`   : plotting primitive kind (:line, :scatter, :hist, :heatmap, ...)
- `xdata`  : x values \\[dimensionless or scaled\\]
- `ydata`  : y values \\[dimensionless or scaled\\]
- `zdata`  : z values if applicable, otherwise `nothing`
- `label`  : legend entry for this series, or `nothing` for no legend
"""
struct SeriesSpec
	kind  :: Symbol
	xdata :: Union{Nothing, AbstractVector{<:Number}}
	ydata :: Union{Nothing, AbstractArray{<:Number}}
	zdata :: Union{Nothing, AbstractArray{<:Number}}
	label :: Union{Nothing, String}
end

"""
	ViewSpec

One plot view / axis system.

All SeriesSpec inside a ViewSpec share the same x/y(/z) axes. The `key`
field encodes the semantic facet this area represents (e.g. \\(i,j\\),
quantity, frequency).

Fields:
- `xaxis`  : x AxisSpec or `nothing` if unused
- `yaxis`  : y AxisSpec or `nothing` if unused
- `zaxis`  : z AxisSpec or `nothing` if unused
- `title`  : view title
- `series` : vector of SeriesSpec
- `key`    : NamedTuple identifying the facet, or empty `NamedTuple` if none
"""
struct ViewSpec
	xaxis  :: Union{Nothing, AxisSpec}
	yaxis  :: Union{Nothing, AxisSpec}
	zaxis  :: Union{Nothing, AxisSpec}
	title  :: String
	series :: Vector{SeriesSpec}
	key    :: NamedTuple
end

"""
	PageSpec

One logical figure / window.

Fields:
- `title`      : optional figure-level title (may be empty)
- `size`       : (width, height) in pixels
- `layout`     : layout spec (:windows, :grid, :tabbed, ...)
- `views`      : vector of ViewSpec values contained in this figure
- `kwargs`     : figure-level backend options (e.g. figsize)
"""
struct PageSpec
	title  :: String
	size   :: Tuple{Int, Int}
	layout :: Symbol
	views  :: Vector{ViewSpec}
	kwargs :: NamedTuple
end

"""
	RenderSpec{S}

Final product of the grammar pipeline for spec type `S`.

Carries only:
- the spec type (the grammar definition),
- a vector of PageSpec payloads.

The rendering backend (Makie) should only see RenderSpec values and must
never touch domain objects or grammar logic.
"""
struct RenderSpec{S <: AbstractPlotSpec}
	spec    :: Type{S}
	figures :: Vector{PageSpec}
end
